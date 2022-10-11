// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#ifndef _SPTAG_SPANN_INDEX_H_
#define _SPTAG_SPANN_INDEX_H_

#include "../Common.h"
#include "../VectorIndex.h"

#include "../Common/CommonUtils.h"
#include "../Common/DistanceUtils.h"
#include "../Common/QueryResultSet.h"
#include "../Common/BKTree.h"
#include "../Common/WorkSpacePool.h"
#include "../Common/FineGrainedLock.h"

#include "../Common/VersionLabel.h"
#include "inc/Helper/SimpleIniReader.h"
#include "inc/Helper/StringConvert.h"
#include "inc/Helper/ThreadPool.h"
#include "inc/Helper/ConcurrentSet.h"
#include "inc/Helper/VectorSetReader.h"

#include "IExtraSearcher.h"
#include "Options.h"
#include "PersistentBuffer.h"

#include <functional>
#include <shared_mutex>
#include <utility>
#include <random>
#include <tbb/concurrent_hash_map.h>

namespace SPTAG
{

    namespace Helper
    {
        class IniReader;
    }

    namespace SPANN
    {
        template<typename T>
    class Index : public VectorIndex
    {
            class AppendAsyncJob : public Helper::ThreadPool::Job
            {
            private:
                VectorIndex* m_index;
                SizeType headID;
                int appendNum;
                std::shared_ptr<std::string> appendPosting;
                std::function<void()> m_callback;
            public:
                AppendAsyncJob(VectorIndex* m_index, SizeType headID, int appendNum, std::shared_ptr<std::string> appendPosting, std::function<void()> p_callback)
                        : m_index(m_index), headID(headID), appendNum(appendNum), appendPosting(std::move(appendPosting)), m_callback(std::move(p_callback)) {}

                ~AppendAsyncJob() {}

                inline void exec(IAbortOperation* p_abort) override {
                    m_index->Append(headID, appendNum, *appendPosting);
                    if (m_callback != nullptr) {
                        m_callback();
                    }
                }
            };

            class ReassignAsyncJob : public SPTAG::Helper::ThreadPool::Job
            {
            private:
                VectorIndex* m_index;
                std::shared_ptr<std::string> vectorContain;
                SizeType VID;
                SizeType HeadPrev;
                uint8_t version;
                std::function<void()> m_callback;
            public:
                ReassignAsyncJob(VectorIndex* m_index,
                                 std::shared_ptr<std::string> vectorContain, SizeType VID, SizeType HeadPrev, uint8_t version, std::function<void()> p_callback)
                        : m_index(m_index),
                          vectorContain(std::move(vectorContain)), VID(VID), HeadPrev(HeadPrev), version(version), m_callback(std::move(p_callback)) {}

                ~ReassignAsyncJob() {}

                void exec(IAbortOperation* p_abort) override {
                    m_index->ProcessAsyncReassign(vectorContain, VID, HeadPrev, version, std::move(m_callback));
                }
            };

            class ThreadPool : public Helper::ThreadPool 
            {
            private:
                std::atomic_uint32_t currentJobs{0};
            public:
                ThreadPool() : Helper::ThreadPool() {}
                
                ~ThreadPool() {}
                
                void init(int numberOfThreads = 1)
                {
                    m_abort.SetAbort(false);
                    for (int i = 0; i < numberOfThreads; i++)
                    {
                        m_threads.emplace_back([this] {
                            Job *j;
                            while (get(j))
                            {
                                try 
                                {
                                    currentJobs++;
                                    j->exec(&m_abort);
                                    currentJobs--;
                                }
                                catch (std::exception& e) {
                                    LOG(Helper::LogLevel::LL_Error, "ThreadPool: exception in %s %s\n", typeid(*j).name(), e.what());
                                }
                                
                                delete j;
                            }
                        });
                    }
                }

                inline uint32_t runningJobs() { return currentJobs; }

                inline bool     allClear()    { return currentJobs == 0 && jobsize() == 0; }
            };

            class Dispatcher
            {
            private:
                std::thread t;

                std::size_t batch;
                std::atomic_bool running{false};
                std::atomic_uint32_t sentAssignment{0};

                Index* m_index;
                std::shared_ptr<PersistentBuffer> m_persistentBuffer;
                std::shared_ptr<ThreadPool> appendThreadPool;
                std::shared_ptr<ThreadPool> reassignThreadPool;

            public:
                Dispatcher(std::shared_ptr<PersistentBuffer> pb, std::size_t batch, std::shared_ptr<ThreadPool> append, std::shared_ptr<ThreadPool> reassign, Index* m_index)
                        : m_persistentBuffer(pb), batch(batch), appendThreadPool(append), reassignThreadPool(reassign), m_index(m_index) {}

                ~Dispatcher() { running = false; t.join(); }

                void dispatch();

                inline void run() { running = true; t = std::thread(&Dispatcher::dispatch, this); }

                inline void stop() { running = false; }

                inline bool allFinished()
                {
                    return sentAssignment == m_persistentBuffer->GetCurrentAssignmentID()
                           && appendThreadPool->allClear()
                           && reassignThreadPool->allClear();
                }

                inline bool allFinishedExceptReassign()
                {
                    return sentAssignment == m_persistentBuffer->GetCurrentAssignmentID()
                           && appendThreadPool->allClear();
                }

                inline bool reassignFinished()
                {
                    return reassignThreadPool->allClear();
                }
            };

            struct EdgeInsert
            {
                EdgeInsert() : headID(INT64_MAX), fullID(INT64_MAX), distance(INT64_MAX), order(0) {}
                uint64_t headID;
                uint64_t fullID;
                float distance;
                char order;
            };

            struct EdgeCompareInsert
                {
                    bool operator()(const EdgeInsert& a, int b) const
                    {
                        return a.headID < b;
                    };

                    bool operator()(int a, const EdgeInsert& b) const
                    {
                        return a < b.headID;
                    };

                    bool operator()(const EdgeInsert& a, const EdgeInsert& b) const
                    {
                        if (a.headID == b.headID)
                        {
                            if (a.distance == b.distance)
                            {
                                return a.fullID < b.fullID;
                            }

                            return a.distance < b.distance;
                        }

                        return a.headID < b.headID;
                    };
                } g_edgeComparerInsert;

        private:
            std::shared_ptr<VectorIndex> m_index;
            std::shared_ptr<std::uint64_t> m_vectorTranslateMap;
            std::unordered_map<std::string, std::string> m_headParameters;
            //std::unique_ptr<std::shared_timed_mutex[]> m_rwLocks;
            COMMON::FineGrainedRWLock m_rwLocks;

            std::unique_ptr<std::atomic_uint32_t[]> m_postingSizes;
            std::atomic_uint64_t m_vectorNum{0};

            std::shared_ptr<IExtraSearcher> m_extraSearcher;
            std::unique_ptr<COMMON::WorkSpacePool<ExtraWorkSpace>> m_workSpacePool;

            Options m_options;

            float(*m_fComputeDistance)(const T* pX, const T* pY, DimensionType length);
            int m_iBaseSquare;
            
            int m_metaDataSize;

            std::shared_ptr<Dispatcher> m_dispatcher;
            std::shared_ptr<PersistentBuffer> m_persistentBuffer;
            std::shared_ptr<Helper::ThreadPool> m_threadPool;
            std::shared_ptr<ThreadPool> m_appendThreadPool;
            std::shared_ptr<ThreadPool> m_reassignThreadPool;

            COMMON::VersionLabel m_versionMap;
            //COMMON::Labelset m_reassignedID;

            tbb::concurrent_hash_map<SizeType, SizeType> m_reassignMap;

            std::atomic_uint32_t m_headMiss{0};
            uint32_t m_appendTaskNum{0};
            uint32_t m_splitTaskNum{0};
            uint32_t m_splitNum{0};
            uint32_t m_theSameHeadNum{0};
            uint32_t m_reAssignNum{0};
            uint32_t m_garbageNum{0};
            std::vector<SizeType> SimplyCountSplit;
            std::vector<SizeType> TheSameHeadSplit;
            std::mutex m_dataAddLock;

        public:
            Index()
            {
                m_fComputeDistance = COMMON::DistanceCalcSelector<T>(m_options.m_distCalcMethod);
                m_iBaseSquare = (m_options.m_distCalcMethod == DistCalcMethod::Cosine) ? COMMON::Utils::GetBase<T>() * COMMON::Utils::GetBase<T>() : 1;
                m_metaDataSize = sizeof(int) + sizeof(uint8_t) + sizeof(float);
            }

            ~Index() {}

            inline std::shared_ptr<VectorIndex> GetMemoryIndex() { return m_index; }
            inline std::shared_ptr<IExtraSearcher> GetDiskIndex() { return m_extraSearcher; }
            inline Options* GetOptions() { return &m_options; }

            inline SizeType GetNumSamples() const { return m_vectorNum.load(); }
            inline DimensionType GetFeatureDim() const { return m_options.m_dim; }
            inline SizeType GetValueSize() const { return m_options.m_dim * sizeof(T); }

            inline int GetCurrMaxCheck() const { return m_options.m_maxCheck; }
            inline int GetNumThreads() const { return m_options.m_iSSDNumberOfThreads; }
            inline DistCalcMethod GetDistCalcMethod() const { return m_options.m_distCalcMethod; }
            inline IndexAlgoType GetIndexAlgoType() const { return IndexAlgoType::SPANN; }
            inline VectorValueType GetVectorValueType() const { return GetEnumValueType<T>(); }
            
            inline float AccurateDistance(const void* pX, const void* pY) const { 
                if (m_options.m_distCalcMethod == DistCalcMethod::L2) return m_fComputeDistance((const T*)pX, (const T*)pY, m_options.m_dim);

                float xy = m_iBaseSquare - m_fComputeDistance((const T*)pX, (const T*)pY, m_options.m_dim);
                float xx = m_iBaseSquare - m_fComputeDistance((const T*)pX, (const T*)pX, m_options.m_dim);
                float yy = m_iBaseSquare - m_fComputeDistance((const T*)pY, (const T*)pY, m_options.m_dim);
                return 1.0f - xy / (sqrt(xx) * sqrt(yy));
            }
            inline float ComputeDistance(const void* pX, const void* pY) const { return m_fComputeDistance((const T*)pX, (const T*)pY, m_options.m_dim); }
            inline bool ContainSample(const SizeType idx) const { return idx < m_options.m_vectorSize; }

            std::shared_ptr<std::vector<std::uint64_t>> BufferSize() const
            {
                std::shared_ptr<std::vector<std::uint64_t>> buffersize(new std::vector<std::uint64_t>);
                auto headIndexBufferSize = m_index->BufferSize();
                buffersize->insert(buffersize->end(), headIndexBufferSize->begin(), headIndexBufferSize->end());
                buffersize->push_back(sizeof(long long) * m_index->GetNumSamples());
                return std::move(buffersize);
            }

            std::shared_ptr<std::vector<std::string>> GetIndexFiles() const
            {
                std::shared_ptr<std::vector<std::string>> files(new std::vector<std::string>);
                auto headfiles = m_index->GetIndexFiles();
                for (auto file : *headfiles) {
                    files->push_back(m_options.m_headIndexFolder + FolderSep + file);
                }
                files->push_back(m_options.m_headIDFile);
                return std::move(files);
            }

            ErrorCode SaveConfig(std::shared_ptr<Helper::DiskPriorityIO> p_configout);
            ErrorCode SaveIndexData(const std::vector<std::shared_ptr<Helper::DiskPriorityIO>>& p_indexStreams);

            ErrorCode LoadConfig(Helper::IniReader& p_reader);
            ErrorCode LoadIndexData(const std::vector<std::shared_ptr<Helper::DiskPriorityIO>>& p_indexStreams);
            ErrorCode LoadIndexDataFromMemory(const std::vector<ByteArray>& p_indexBlobs);

            ErrorCode BuildIndex(const void* p_data, SizeType p_vectorNum, DimensionType p_dimension, bool p_normalized = false);
            ErrorCode BuildIndex(bool p_normalized = false);
            ErrorCode SearchIndex(QueryResult &p_query, bool p_searchDeleted = false) const;
            ErrorCode DebugSearchDiskIndex(QueryResult& p_query, int p_subInternalResultNum, int p_internalResultNum,
                SearchStats* p_stats = nullptr, std::set<int>* truth = nullptr, std::map<int, std::set<int>>* found = nullptr);
            ErrorCode UpdateIndex();

            ErrorCode SetParameter(const char* p_param, const char* p_value, const char* p_section = nullptr);
            std::string GetParameter(const char* p_param, const char* p_section = nullptr) const;

            inline const void* GetSample(const SizeType idx) const { return nullptr; }
            inline SizeType GetNumDeleted() const { return 0; }
            inline bool NeedRefine() const { return false; }
            inline bool CheckIdDeleted(const SizeType& p_id) { return m_versionMap.Contains(p_id); }
            inline bool CheckVersionValid(const SizeType& p_id, const uint8_t version) {return m_versionMap.GetVersion(p_id) == version;}

            ErrorCode RefineSearchIndex(QueryResult &p_query, bool p_searchDeleted = false) const { return ErrorCode::Undefined; }
            ErrorCode SearchTree(QueryResult& p_query) const { return ErrorCode::Undefined; }
            ErrorCode AddIndex(const void* p_data, SizeType p_vectorNum, DimensionType p_dimension, std::shared_ptr<MetadataSet> p_metadataSet, bool p_withMetaIndex = false, bool p_normalized = false);
            ErrorCode AddIndexId(const void* p_data, SizeType p_vectorNum, DimensionType p_dimension, int& beginHead, int& endHead)  { return ErrorCode::Undefined; }
            ErrorCode AddIndexIdx(SizeType begin, SizeType end) { return ErrorCode::Undefined; }
            ErrorCode DeleteIndex(const void* p_vectors, SizeType p_vectorNum) { return ErrorCode::Undefined; }
            ErrorCode DeleteIndex(const SizeType& p_id);
            ErrorCode RefineIndex(const std::vector<std::shared_ptr<Helper::DiskPriorityIO>>& p_indexStreams, IAbortOperation* p_abort) { return ErrorCode::Undefined; }
            ErrorCode RefineIndex(std::shared_ptr<VectorIndex>& p_newIndex) { return ErrorCode::Undefined; }
            
        private:
            bool CheckHeadIndexType();
            void SelectHeadAdjustOptions(int p_vectorCount);
            int SelectHeadDynamicallyInternal(const std::shared_ptr<COMMON::BKTree> p_tree, int p_nodeID, const Options& p_opts, std::vector<int>& p_selected);
            void SelectHeadDynamically(const std::shared_ptr<COMMON::BKTree> p_tree, int p_vectorCount, std::vector<int>& p_selected);
            bool SelectHead(std::shared_ptr<Helper::VectorSetReader>& p_reader);

            ErrorCode BuildIndexInternal(std::shared_ptr<Helper::VectorSetReader>& p_reader);

            ErrorCode Append(SizeType headID, int appendNum, std::string& appendPosting);
            ErrorCode Split(const SizeType headID, int appendNum, std::string& appendPosting);
            ErrorCode ReAssign(SizeType headID, std::vector<std::string>& postingLists, std::vector<SizeType>& newHeadsID);
            void ReAssignVectors(std::map<SizeType, T*>& reAssignVectors, std::map<SizeType, SizeType>& HeadPrevs, std::map<SizeType, uint8_t>& versions);
            void ReAssignUpdate(const std::shared_ptr<std::string>&, SizeType VID, SizeType HeadPrev, uint8_t version);

        public:
            inline void AppendAsync(SizeType headID, int appendNum, std::shared_ptr<std::string> appendPosting, std::function<void()> p_callback=nullptr)
            {
                auto* curJob = new AppendAsyncJob(this, headID, appendNum, std::move(appendPosting), p_callback);
                m_appendThreadPool->add(curJob);
            }

            inline void ReassignAsync(std::shared_ptr<std::string> vectorContain, SizeType VID, SizeType HeadPrev, uint8_t version, std::function<void()> p_callback=nullptr)
            {   
                auto* curJob = new ReassignAsyncJob(this, std::move(vectorContain), VID, HeadPrev, version, p_callback);
                m_reassignThreadPool->add(curJob);
            }

            void ProcessAsyncReassign(std::shared_ptr<std::string> vectorContain, SizeType VID, SizeType HeadPrev, uint8_t version, std::function<void()> p_callback);

            bool AllFinished() {return m_dispatcher->allFinished();}

            bool AllFinishedExceptReassign() {return m_dispatcher->allFinishedExceptReassign();}

            bool ReassignFinished() {return m_dispatcher->reassignFinished();}

            void ForceCompaction() {if (m_options.m_useKV) m_extraSearcher->ForceCompaction();}

            int getSplitTimes() {return m_splitNum;}

            int getHeadMiss() {return m_headMiss.load();}

            int getSameHead() {return m_theSameHeadNum;}

            int getReassignNum() {return m_reAssignNum;}
            
            int getGarbageNum() {return m_garbageNum;}

            void printSplitStatus() 
            {
                for (int i = 0; i < 13; i++)
                {
                    LOG(Helper::LogLevel::LL_Info, "%d ~ %d\t", i*10, (i+1)*10-1);
                }
                LOG(Helper::LogLevel::LL_Info, "\n");
                for (int i = 0; i < 13; i++)
                {
                    LOG(Helper::LogLevel::LL_Info, "%6d\t", SimplyCountSplit[i]);
                }
                LOG(Helper::LogLevel::LL_Info, "\n");
                return;
            }

            void printSplitTheSameHeadStatus()
            {
                int postingNum = m_index->GetNumSamples();
                std::map<int, int> count;
                for (int i = 0; i < postingNum; i++)
                {
                    if(count.find(TheSameHeadSplit[i]) == count.end()) {
                        count[TheSameHeadSplit[i]] = 0;
                    }
                    count[TheSameHeadSplit[i]]++;
                }
                for (auto it = count.begin(); it != count.end(); it++)
                {
                    LOG(Helper::LogLevel::LL_Info, "TheSameHeadSplit %d times: %d nodes\n", it->first, it->second);
                }
            }

            void UpdateStop()
            {
                m_persistentBuffer->StopPB();
                m_dispatcher->stop();
            }

            void Rebuild(std::shared_ptr<Helper::VectorSetReader>& p_reader, SizeType upperBound = -1)
            {
                auto fullVectors = p_reader->GetVectorSet();
                int curCount;
                if (upperBound == -1) {
                    curCount = fullVectors->Count();
                } else {
                    curCount = upperBound;
                }
                LOG(Helper::LogLevel::LL_Info, "Rebuild SSD Index.\n");
                std::vector<EdgeInsert> selections(static_cast<size_t>(curCount)* m_options.m_replicaCount);

                std::vector<int> replicaCount(curCount, 0);
                std::vector<std::atomic_int> postingListSize(m_index->GetNumSamples());
                for (auto& pls : postingListSize) pls = 0;
                LOG(Helper::LogLevel::LL_Info, "Preparation done, start candidate searching.\n");

                std::vector<std::thread> threads;
                threads.reserve(64);

                std::atomic_int nextFullID(0);
                std::atomic_size_t rngFailedCountTotal(0);

                for (int tid = 0; tid < 64; ++tid)
                {
                    threads.emplace_back([&, tid]()
                        {
                            COMMON::QueryResultSet<T> resultSet(NULL, m_options.m_internalResultNum);

                            size_t rngFailedCount = 0;

                            while (true)
                            {
                                int fullID = nextFullID.fetch_add(1);
                                if (fullID >= curCount)
                                {
                                    break;
                                }

                                T* buffer = reinterpret_cast<T*>(fullVectors->GetVector(fullID));
                                resultSet.SetTarget(buffer);
                                resultSet.Reset();

                                m_index->SearchIndex(resultSet);

                                size_t selectionOffset = static_cast<size_t>(fullID)* m_options.m_replicaCount;

                                BasicResult* queryResults = resultSet.GetResults();
                                for (int i = 0; i < m_options.m_internalResultNum && replicaCount[fullID] < m_options.m_replicaCount; ++i)
                                {
                                    if (queryResults[i].VID == -1)
                                    {
                                        break;
                                    }

                                    // RNG Check.
                                    bool rngAccpeted = true;
                                    for (int j = 0; j < replicaCount[fullID]; ++j)
                                    {
                                        // VQANNSearch::QueryResultSet<ValueType> resultSet(NULL, candidateNum);

                                        float nnDist = m_index->ComputeDistance(
                                            m_index->GetSample(queryResults[i].VID),
                                            m_index->GetSample(selections[selectionOffset + j].headID));

                                        // LOG(Helper::LogLevel::LL_Info,  "NNDist: %f Original: %f\n", nnDist, queryResults[i].Score);
                                        if (nnDist <= queryResults[i].Dist)
                                        {
                                            rngAccpeted = false;
                                            break;
                                        }
                                    }

                                    if (!rngAccpeted)
                                    {
                                        ++rngFailedCount;
                                        continue;
                                    }

                                    ++postingListSize[queryResults[i].VID];
                                    selections[selectionOffset + replicaCount[fullID]].headID = queryResults[i].VID;
                                    selections[selectionOffset + replicaCount[fullID]].fullID = fullID;
                                    selections[selectionOffset + replicaCount[fullID]].distance = queryResults[i].Dist;
                                    selections[selectionOffset + replicaCount[fullID]].order = (char)replicaCount[fullID];
                                    ++replicaCount[fullID];
                                }
                            }

                            rngFailedCountTotal += rngFailedCount;
                        });
                }

                for (int tid = 0; tid < 64; ++tid)
                {
                    threads[tid].join();
                }

                LOG(Helper::LogLevel::LL_Info, "Searching replicas ended. RNG failed count: %llu\n", static_cast<uint64_t>(rngFailedCountTotal.load()));

                std::sort(selections.begin(), selections.end(), g_edgeComparerInsert);

                int postingSizeLimit = INT_MAX;
                if (m_options.m_postingPageLimit > 0)
                {
                    postingSizeLimit = static_cast<int>(m_options.m_postingPageLimit * 4096/ (fullVectors->PerVectorDataSize() + sizeof(int)));
                }

                LOG(Helper::LogLevel::LL_Info, "Posting size limit: %d\n", postingSizeLimit);

                {
                    std::vector<int> replicaCountDist(m_options.m_replicaCount + 1, 0);
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        ++replicaCountDist[replicaCount[i]];
                    }

                    LOG(Helper::LogLevel::LL_Info, "Before Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        LOG(Helper::LogLevel::LL_Info, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
                    }
                }
                for (int i = 0; i < postingListSize.size(); ++i)
                {
                    std::size_t selectIdx = std::lower_bound(selections.begin(), selections.end(), i, g_edgeComparerInsert) - selections.begin();
                    if (postingListSize[i] <= postingSizeLimit) {
                        continue;
                    }
                    for (size_t dropID = postingSizeLimit; dropID < postingListSize[i]; ++dropID) {
                        int fullID = selections[selectIdx + dropID].fullID;
                        --replicaCount[fullID];
                    }

                    postingListSize[i] = postingSizeLimit;
                }
                std::string postinglist;
                for (int id = 0; id < postingListSize.size(); id++) 
                {
                    postinglist.resize(0);
                    postinglist.clear();
                    std::size_t selectIdx = std::lower_bound(selections.begin(), selections.end(), id, g_edgeComparerInsert)
                                            - selections.begin();
                    for (int j = 0; j < postingListSize[id]; ++j) {
                        if (selections[selectIdx].headID != id) {
                            LOG(Helper::LogLevel::LL_Error, "Selection ID NOT MATCH\n");
                            exit(1);
                        }
                        float distance = selections[selectIdx].distance;
                        int fullID = selections[selectIdx++].fullID;
                        uint8_t version = 0;
                        m_versionMap.UpdateVersion(fullID, 0);
                        size_t dim = fullVectors->Dimension();
                        // First Vector ID, then Vector
                        postinglist += Helper::Convert::Serialize<int>(&fullID, 1);
                        postinglist += Helper::Convert::Serialize<uint8_t>(&version, 1);
                        postinglist += Helper::Convert::Serialize<float>(&distance, 1);
                        postinglist += Helper::Convert::Serialize<T>(fullVectors->GetVector(fullID), dim);
                    }
                    m_extraSearcher->AddIndex(id, postinglist);
                }
                for(int id = 0; id < postingListSize.size(); id++)
                {
                    m_postingSizes[id].store(postingListSize[id]);
                }
            }

            void QuantifyAssumptionBrokenTotally()
            {
                std::atomic_int assumptionBrokenNum = 0;
                std::atomic_int assumption64BrokenNum = 0;
                std::atomic_int deleted = 0;
                int vectorNum = m_vectorNum.load();
                std::vector<std::set<SizeType>> vectorHeadMap(vectorNum);
                std::vector<bool> vectorFoundMap(vectorNum);
                std::vector<std::string> vectorIdValueMap(vectorNum);
                for (int i = 0; i < m_index->GetNumSamples(); i++) {
                    std::string postingList;
                    if (!m_index->ContainSample(i)) continue;
                    m_extraSearcher->SearchIndex(i, postingList);
                    int postVectorNum = postingList.size() / (m_options.m_dim * sizeof(T) + m_metaDataSize);
                    uint8_t* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                    for (int j = 0; j < postVectorNum; j++) {
                        uint8_t* vectorId = postingP + j * (m_options.m_dim * sizeof(T) + m_metaDataSize);
                        SizeType vid = *(reinterpret_cast<SizeType*>(vectorId));
                        if (m_versionMap.Contains(vid)) continue;
                        vectorHeadMap[vid].insert(i);
                        if (vectorFoundMap[vid]) continue;
                        vectorFoundMap[vid] = true;
                        vectorIdValueMap[vid] = Helper::Convert::Serialize<uint8_t>(vectorId + m_metaDataSize, m_options.m_dim * sizeof(T));
                    }
                }
                LOG(Helper::LogLevel::LL_Info, "Begin Searching\n");
                #pragma omp parallel for num_threads(32)
                for (int vid = 0; vid < vectorNum; vid++) {
                    if (m_versionMap.Contains(vid)) {
                        deleted++;
                        continue;
                    }
                    COMMON::QueryResultSet<T> headCandidates(NULL, 64);
                    headCandidates.SetTarget(reinterpret_cast<T*>(&vectorIdValueMap[vid].front()));
                    headCandidates.Reset();
                    m_index->SearchIndex(headCandidates);
                    int replicaCount = 0;
                    BasicResult* queryResults = headCandidates.GetResults();
                    std::vector<EdgeInsert> selections(static_cast<size_t>(m_options.m_replicaCount));
                    std::set<SizeType>HeadMap;
                    for (int i = 0; i < headCandidates.GetResultNum() && replicaCount < m_options.m_replicaCount; ++i) {
                        if (queryResults[i].VID == -1) {
                            break;
                        }
                        HeadMap.insert(queryResults[i].VID);
                        // RNG Check.
                        bool rngAccpeted = true;
                        for (int j = 0; j < replicaCount; ++j) {
                            float nnDist = m_index->ComputeDistance(
                                                        m_index->GetSample(queryResults[i].VID),
                                                        m_index->GetSample(selections[j].headID));
                            if (nnDist <= queryResults[i].Dist) {
                                rngAccpeted = false;
                                break;
                            }
                        }
                        if (!rngAccpeted)
                            continue;

                        selections[replicaCount].headID = queryResults[i].VID;
                        ++replicaCount;
                    }
                    for (int i = 0; i < replicaCount; i++) {
                        if (!vectorHeadMap[vid].count(selections[i].headID)) {
                            assumptionBrokenNum++;
                            break;
                        }
                    }
                    for (auto vectorHead : vectorHeadMap[vid]) {
                        if (!HeadMap.count(vectorHead)) {
                            assumption64BrokenNum++;
                            break;
                        }
                    }
                }
                LOG(Helper::LogLevel::LL_Info, "After Split %d times, %d total vectors, %d assumption broken vectors, %d vectors are not in 64 nearest postings\n", m_splitNum, vectorNum - deleted.load(), assumptionBrokenNum.load(), assumption64BrokenNum.load());
            }

            //headCandidates: search data structrue for "vid" vector
            //headID: the head vector that stands for vid
            bool IsAssumptionBroken(SizeType headID, COMMON::QueryResultSet<T>& headCandidates, SizeType vid)
            {
                m_index->SearchIndex(headCandidates);
                int replicaCount = 0;
                BasicResult* queryResults = headCandidates.GetResults();
                std::vector<EdgeInsert> selections(static_cast<size_t>(m_options.m_replicaCount));
                for (int i = 0; i < headCandidates.GetResultNum() && replicaCount < m_options.m_replicaCount; ++i) {
                    if (queryResults[i].VID == -1) {
                        break;
                    }
                    // RNG Check.
                    bool rngAccpeted = true;
                    for (int j = 0; j < replicaCount; ++j) {
                        float nnDist = m_index->ComputeDistance(
                                                    m_index->GetSample(queryResults[i].VID),
                                                    m_index->GetSample(selections[j].headID));
                        if (nnDist <= queryResults[i].Dist) {
                            rngAccpeted = false;
                            break;
                        }
                    }
                    if (!rngAccpeted)
                        continue;

                    selections[replicaCount].headID = queryResults[i].VID;
                    // LOG(Helper::LogLevel::LL_Info, "head:%d\n", queryResults[i].VID);
                    if (selections[replicaCount].headID == headID) return false;
                    ++replicaCount;
                }
                return true;
            }

            //Measure that in "headID" posting list, how many vectors break their assumption
            int QuantifyAssumptionBroken(SizeType headID, std::string& postingList, SizeType SplitHead, std::vector<SizeType>& newHeads, std::set<int>& brokenID, int topK = 0, float ratio = 1.0)
            {
                int assumptionBrokenNum = 0;
                int m_vectorInfoSize = sizeof(T) * m_options.m_dim + m_metaDataSize;
                int postVectorNum = postingList.size() / m_vectorInfoSize;
                uint8_t* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                float minDist;
                float maxDist;
                float avgDist = 0;
                std::vector<float> distanceSet;
                //#pragma omp parallel for num_threads(32)
                for (int j = 0; j < postVectorNum; j++) {
                    uint8_t* vectorId = postingP + j * m_vectorInfoSize;
                    SizeType vid = *(reinterpret_cast<int*>(vectorId));
                    uint8_t version = *(reinterpret_cast<uint8_t*>(vectorId + sizeof(int)));
                    float_t dist = *(reinterpret_cast<float*>(vectorId + sizeof(int) + sizeof(uint8_t)));
                    // if (dist < Epsilon) LOG(Helper::LogLevel::LL_Info, "head found: vid: %d, head: %d\n", vid, headID);
                    avgDist += dist;
                    distanceSet.push_back(dist);
                    if (CheckIdDeleted(vid) || !CheckVersionValid(vid, version)) continue;
                    COMMON::QueryResultSet<T> headCandidates(NULL, 64);
                    headCandidates.SetTarget(reinterpret_cast<T*>(vectorId + m_metaDataSize));
                    headCandidates.Reset();
                    if (brokenID.find(vid) == brokenID.end() && IsAssumptionBroken(headID, headCandidates, vid)) {
                        /*
                        float_t headDist = m_index->ComputeDistance(headCandidates.GetTarget(), m_index->GetSample(SplitHead));
                        float_t newHeadDist_1 = m_index->ComputeDistance(headCandidates.GetTarget(), m_index->GetSample(newHeads[0]));
                        float_t newHeadDist_2 = m_index->ComputeDistance(headCandidates.GetTarget(), m_index->GetSample(newHeads[1]));

                        float_t splitDist = m_index->ComputeDistance(m_index->GetSample(SplitHead), m_index->GetSample(headID));

                        float_t headToNewHeadDist_1 = m_index->ComputeDistance(m_index->GetSample(headID), m_index->GetSample(newHeads[0]));
                        float_t headToNewHeadDist_2 = m_index->ComputeDistance(m_index->GetSample(headID), m_index->GetSample(newHeads[1]));

                        LOG(Helper::LogLevel::LL_Info, "broken vid to head distance: %f, to split head distance: %f\n", dist, headDist);
                        LOG(Helper::LogLevel::LL_Info, "broken vid to new head 1 distance: %f, to new head 2 distance: %f\n", newHeadDist_1, newHeadDist_2);
                        LOG(Helper::LogLevel::LL_Info, "head to spilit head distance: %f\n", splitDist);
                        LOG(Helper::LogLevel::LL_Info, "head to new head 1 distance: %f, to new head 2 distance: %f\n", headToNewHeadDist_1, headToNewHeadDist_2);
                        */
                        assumptionBrokenNum++;
                        brokenID.insert(vid);
                    }
                }
                
                if (assumptionBrokenNum != 0) {
                    std::sort(distanceSet.begin(), distanceSet.end());
                    minDist = distanceSet[1];
                    maxDist = distanceSet[distanceSet.size() - 1];
                    // LOG(Helper::LogLevel::LL_Info, "distance: min: %f, max: %f, avg: %f, 50th: %f\n", minDist, maxDist, avgDist/postVectorNum, distanceSet[distanceSet.size() * 0.5]);
                    LOG(Helper::LogLevel::LL_Info, "assumption broken num: %d\n", assumptionBrokenNum);
                    float_t splitDist = m_index->ComputeDistance(m_index->GetSample(SplitHead), m_index->GetSample(headID));

                    float_t headToNewHeadDist_1 = m_index->ComputeDistance(m_index->GetSample(headID), m_index->GetSample(newHeads[0]));
                    float_t headToNewHeadDist_2 = m_index->ComputeDistance(m_index->GetSample(headID), m_index->GetSample(newHeads[1]));

                    LOG(Helper::LogLevel::LL_Info, "head to spilt head distance: %f/%d/%.2f\n", splitDist, topK, ratio);
                    LOG(Helper::LogLevel::LL_Info, "head to new head 1 distance: %f, to new head 2 distance: %f\n", headToNewHeadDist_1, headToNewHeadDist_2);
                }
                
                return assumptionBrokenNum;
            }

            int QuantifySplitCaseA(std::vector<SizeType>& newHeads, std::vector<std::string>& postingLists, SizeType SplitHead, int split_order, std::set<int>& brokenID)
            {
                int assumptionBrokenNum = 0;
                assumptionBrokenNum += QuantifyAssumptionBroken(newHeads[0], postingLists[0], SplitHead, newHeads, brokenID);
                assumptionBrokenNum += QuantifyAssumptionBroken(newHeads[1], postingLists[1], SplitHead, newHeads, brokenID);
                int vectorNum = (postingLists[0].size() + postingLists[1].size()) / (sizeof(T) * m_options.m_dim + m_metaDataSize);
                // LOG(Helper::LogLevel::LL_Info, "After Split%d, Top0 nearby posting lists, caseA : %d/%d\n", split_order, assumptionBrokenNum, vectorNum);
                return assumptionBrokenNum;
            }

            //Measure that around "headID", how many vectors break their assumption
            //"headID" is the head vector before split
            void QuantifySplitCaseB(SizeType headID, std::vector<SizeType>& newHeads, SizeType SplitHead, int split_order, int assumptionBrokenNum_top0, std::set<int>& brokenID)
            {
                auto headVector = reinterpret_cast<const T*>(m_index->GetSample(headID));
                COMMON::QueryResultSet<T> nearbyHeads(NULL, 128);
                nearbyHeads.SetTarget(headVector);
                nearbyHeads.Reset();
                std::vector<std::string> postingLists;
                m_index->SearchIndex(nearbyHeads);
                std::string postingList;
                BasicResult* queryResults = nearbyHeads.GetResults();
                int topk = 8;
                int assumptionBrokenNum = assumptionBrokenNum_top0;
                int assumptionBrokenNum_topK = assumptionBrokenNum_top0;
                int i;
                int containedHead = 0;
                if (assumptionBrokenNum_top0 != 0) containedHead++;
                int vectorNum = 0;
                float furthestDist = 0;
                for (i = 0; i < nearbyHeads.GetResultNum(); i++) {
                    if (queryResults[i].VID == -1) {
                        break;
                    }
                    furthestDist = queryResults[i].Dist;
                    if (i == topk) {
                        LOG(Helper::LogLevel::LL_Info, "After Split%d, Top%d nearby posting lists, caseB : %d in %d/%d\n", split_order, i, assumptionBrokenNum, containedHead, vectorNum);
                        LOG(Helper::LogLevel::LL_Info, "After Split%d, Top%d Dist: %f\n", split_order, i, furthestDist);
                        topk *= 2;
                    }
                    if (queryResults[i].VID == newHeads[0] || queryResults[i].VID == newHeads[1]) continue;
                    m_extraSearcher->SearchIndex(queryResults[i].VID, postingList);
                    vectorNum += postingList.size() / (sizeof(T) * m_options.m_dim + m_metaDataSize);
                    int tempNum = QuantifyAssumptionBroken(queryResults[i].VID, postingList, SplitHead, newHeads, brokenID, i, queryResults[i].Dist/queryResults[1].Dist);
                    assumptionBrokenNum += tempNum;
                    if (tempNum != 0) containedHead++;
                }
                LOG(Helper::LogLevel::LL_Info, "After Split%d, Top%d nearby posting lists, caseB : %d in %d/%d\n", split_order, i, assumptionBrokenNum, containedHead, vectorNum);
                LOG(Helper::LogLevel::LL_Info, "After Split%d, Top%d Dist: %f\n", split_order, i, furthestDist);
            }

            void QuantifySplit(SizeType headID, std::vector<std::string>& postingLists, std::vector<SizeType>& newHeads, SizeType SplitHead, int split_order)
            {
                std::set<int> brokenID;
                brokenID.clear();
                // LOG(Helper::LogLevel::LL_Info, "Split Quantify: %d, head1:%d, head2:%d\n", split_order, newHeads[0], newHeads[1]);
                int assumptionBrokenNum = QuantifySplitCaseA(newHeads, postingLists, SplitHead, split_order, brokenID);
                QuantifySplitCaseB(headID, newHeads, SplitHead, split_order, assumptionBrokenNum, brokenID);
            }

            bool CheckIsInMiddleOfTwoHead(T* data, std::vector<SizeType>& newHeads, SizeType currentHead, float_t currentHeadDist, SizeType splitHead)
            {
                float_t splitHeadDist = m_index->ComputeDistance(data, m_index->GetSample(splitHead));

                float_t newHeadToHeadDist_1 = m_index->ComputeDistance(m_index->GetSample(newHeads[0]), m_index->GetSample(currentHead));
                float_t newHeadToVectorDist_1 = m_index->ComputeDistance(m_index->GetSample(newHeads[0]), data);

                float_t newHeadToHeadDist_2 = m_index->ComputeDistance(m_index->GetSample(newHeads[1]), m_index->GetSample(currentHead));
                float_t newHeadToVectorDist_2 = m_index->ComputeDistance(m_index->GetSample(newHeads[1]), data);
                // if (newHeadToHeadDist_1 >= newHeadToVectorDist_1 || newHeadToHeadDist_2 >= newHeadToVectorDist_2) return true;

                if (newHeadToVectorDist_1 < currentHeadDist || newHeadToVectorDist_2 < currentHeadDist) return true;

                // LOG(Helper::LogLevel::LL_Info, "NewHeadToHeadDist_1: %f, NewHeadToHeadDist_2: %f, newHeadToVectorDist_1: %f, newHeadToVectorDist_2: %f,  currentHeadDist: %f, splitHeadDist: %f\n", 
                //  newHeadToHeadDist_1, newHeadToHeadDist_2, newHeadToVectorDist_1, newHeadToVectorDist_2, currentHeadDist, splitHeadDist);
                return false;
            }

            bool CheckIsReallyNeedReassign(T* data, std::vector<SizeType>& newHeads, SizeType currentHead)
            {
                COMMON::QueryResultSet<T> p_queryResults(NULL, m_options.m_internalResultNum);
                p_queryResults.SetTarget(reinterpret_cast<T*>(data));
                p_queryResults.Reset();
                m_index->SearchIndex(p_queryResults);
                BasicResult* queryResults = p_queryResults.GetResults();
                for (int i = 0; i < p_queryResults.GetResultNum(); ++i) {
                    if (queryResults[i].VID == -1 || queryResults[i].VID == currentHead) {
                        break;
                    }
                    for (auto headID : newHeads) {
                        if (queryResults[i].VID == headID) return true;
                    }
                }
                return false;
            }

            bool CheckIsNeedReassign(std::vector<SizeType>& newHeads, T* data, SizeType splitHead, float_t headToSplitHeadDist, float_t currentHeadDist, bool isInSplitHead, SizeType currentHead)
            {
                
                float_t headDist = m_index->ComputeDistance(data, m_index->GetSample(splitHead));
                float_t newHeadDist_1 = m_index->ComputeDistance(data, m_index->GetSample(newHeads[0]));
                float_t newHeadDist_2 = m_index->ComputeDistance(data, m_index->GetSample(newHeads[1]));

                // float_t newHeadToHeadDist_1 = m_index->ComputeDistance(m_index->GetSample(newHeads[0]), m_index->GetSample(currentHead));
                // float_t newHeadToHeadDist_2 = m_index->ComputeDistance(m_index->GetSample(newHeads[1]), m_index->GetSample(currentHead));

                // float_t newHeadToSplitHeadDist_1 = m_index->ComputeDistance(m_index->GetSample(newHeads[0]), m_index->GetSample(splitHead));
                // float_t newHeadToSplitHeadDist_2 = m_index->ComputeDistance(m_index->GetSample(newHeads[1]), m_index->GetSample(splitHead));
                // float_t newHeadToNewHeadDist = m_index->ComputeDistance(m_index->GetSample(newHeads[0]), m_index->GetSample(newHeads[1]));


                // if (headDist * headDist >= (currentHeadDist * currentHeadDist + headToSplitHeadDist * headToSplitHeadDist) ) return false;

                if (isInSplitHead) {
                    if (headDist >= currentHeadDist) return false;
                } else {
                    if (headDist <= newHeadDist_1 && headDist <= newHeadDist_2) return false;
                    if (currentHeadDist <= newHeadDist_1 && currentHeadDist <= newHeadDist_2) return false;
                }

                // LOG(Helper::LogLevel::LL_Info, "newHeadDist_1: %f, newHeadDist_2: %f, currentHeadDist: %f, toSplitHeadDist: %f, headToSplitHeadDist: %f, headToNewHeadDist_1: %f, headToNewHeadDist_2: %f newHeadToSplitHeadDist_1: %f, newHeadToSplitHeadDist_2: %f, newHeadToNewHeadDist: %f\n",
                //     newHeadDist_1, newHeadDist_2, currentHeadDist, headDist, headToSplitHeadDist, newHeadToHeadDist_1, newHeadToHeadDist_2,
                //     newHeadToSplitHeadDist_1, newHeadToSplitHeadDist_2, newHeadToNewHeadDist);
                // if ((headToSplitHeadDist + newHeadToSplitHeadDist_1 ) > newHeadToHeadDist_1) LOG(Helper::LogLevel::LL_Info, "Head_1: True\n");
                // else LOG(Helper::LogLevel::LL_Info, "Head_1: False\n");
                // if ((headToSplitHeadDist + newHeadToSplitHeadDist_2 ) > newHeadToHeadDist_2) LOG(Helper::LogLevel::LL_Info, "Head_2: True\n");
                // else LOG(Helper::LogLevel::LL_Info, "Head_2: False\n");
                // COMMON::QueryResultSet<T> headCandidates(NULL, 64);
                // headCandidates.SetTarget(reinterpret_cast<T*>(data));
                // headCandidates.Reset();
                // if (IsAssumptionBroken(currentHead, headCandidates, 0)) LOG(Helper::LogLevel::LL_Info, "Really Need Reassign\n");
                // if (currentHeadDist > headDist)
                // {
                //     std::string postingList;
                //     m_extraSearcher->SearchIndex(newHeads[0], postingList);
                //     int vectorInfoSize = m_options.m_dim * sizeof(T) + m_metaDataSize;
                //     size_t postVectorNum = postingList.size() / vectorInfoSize;
                //     auto* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                //     for (int j = 0; j < postVectorNum; j++) {
                //         uint8_t* vectorId = postingP + j * vectorInfoSize;
                //         SizeType vid = *(reinterpret_cast<SizeType*>(vectorId));
                //         float dist = *(reinterpret_cast<float*>(vectorId + sizeof(int) + sizeof(uint8_t)));
                //         if (dist == headDist) {
                //             LOG(Helper::LogLevel::LL_Info, "It is in splithead\n");
                //             return true;
                //         }
                //     }
                //     m_extraSearcher->SearchIndex(newHeads[1], postingList);
                //     postVectorNum = postingList.size() / vectorInfoSize;
                //     postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                //     for (int j = 0; j < postVectorNum; j++) {
                //         uint8_t* vectorId = postingP + j * vectorInfoSize;
                //         SizeType vid = *(reinterpret_cast<SizeType*>(vectorId));
                //         float dist = *(reinterpret_cast<float*>(vectorId + sizeof(int) + sizeof(uint8_t)));
                //         if (dist == headDist) {
                //             LOG(Helper::LogLevel::LL_Info, "It is in splithead\n");
                //             return true;
                //         }
                //     }
                //     LOG(Helper::LogLevel::LL_Info, "It is not in splithead\n");
                // }

                return true;
            }
        };
    } // namespace SPANN
} // namespace SPTAG

#endif // _SPTAG_SPANN_INDEX_H_
