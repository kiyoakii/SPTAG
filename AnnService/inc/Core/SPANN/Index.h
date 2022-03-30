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

#include "../Common/Labelset.h"
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
                std::vector<SizeType>& newHeads;
                bool check;
                SizeType oldVID;
                std::function<void()> m_callback;
            public:
                ReassignAsyncJob(VectorIndex* m_index,
                                 std::shared_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check,
                                 SizeType oldVID, std::function<void()> p_callback)
                        : m_index(m_index),
                          vectorContain(std::move(vectorContain)), VID(VID), newHeads(newHeads), check(check), oldVID(oldVID), m_callback(std::move(p_callback)) {}

                ~ReassignAsyncJob() {}

                void exec(IAbortOperation* p_abort) override {
                    m_index->ProcessAsyncReassign(vectorContain, VID, newHeads, check, oldVID, std::move(m_callback));
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
                           && appendThreadPool->runningJobs() == 0
                           && reassignThreadPool->runningJobs() == 0;
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

        private:
            std::shared_ptr<VectorIndex> m_index;
            std::shared_ptr<std::uint64_t> m_vectorTranslateMap;
            std::unordered_map<std::string, std::string> m_headParameters;
            std::unique_ptr<std::shared_timed_mutex[]> m_rwLocks;
            std::unique_ptr<std::atomic_uint32_t[]> m_postingSizes;
            std::atomic_uint64_t m_vectorNum{0};

            std::shared_ptr<IExtraSearcher> m_extraSearcher;
            std::unique_ptr<COMMON::WorkSpacePool<ExtraWorkSpace>> m_workSpacePool;

            Options m_options;

            float(*m_fComputeDistance)(const T* pX, const T* pY, DimensionType length);
            int m_iBaseSquare;

            std::shared_ptr<Dispatcher> m_dispatcher;
            std::shared_ptr<PersistentBuffer> m_persistentBuffer;
            std::shared_ptr<Helper::ThreadPool> m_threadPool;
            std::shared_ptr<ThreadPool> m_appendThreadPool;
            std::shared_ptr<ThreadPool> m_reassignThreadPool;

            COMMON::Labelset m_deletedID;
            COMMON::Labelset m_reassignedID;
            std::atomic_uint32_t m_headMiss{0};
            uint32_t m_appendTaskNum{0};
            uint32_t m_splitTaskNum{0};
            uint32_t m_splitNum{0};
            std::mutex m_dataAddLock;

        public:
            Index()
            {
                m_fComputeDistance = COMMON::DistanceCalcSelector<T>(m_options.m_distCalcMethod);
                m_iBaseSquare = (m_options.m_distCalcMethod == DistCalcMethod::Cosine) ? COMMON::Utils::GetBase<T>() * COMMON::Utils::GetBase<T>() : 1;

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
            inline bool CheckIdDeleted(const SizeType& p_id) { return m_deletedID.Contains(p_id); }

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

            ErrorCode Append(SizeType headID, int appendNum, std::string& appendPosting, SizeType oldVID);
            ErrorCode Split(const SizeType headID, int appendNum, std::string& appendPosting);
            ErrorCode ReAssign(SizeType headID, std::vector<std::string>& postingLists, std::vector<SizeType>& newHeadsID);
            void ReAssignVectors(std::map<SizeType, T*>& reAssignVectors, std::vector<SizeType>& newHeadsID, bool check=false);
            void ReAssignUpdate(const std::shared_ptr<std::string>&, SizeType VID, std::vector<SizeType>&, bool check = false, SizeType oldVID = 0);

        public:
            inline void AppendAsync(SizeType headID, int appendNum, std::shared_ptr<std::string> appendPosting, std::function<void()> p_callback=nullptr)
            {
                auto* curJob = new AppendAsyncJob(this, headID, appendNum, std::move(appendPosting), p_callback);
                m_appendThreadPool->add(curJob);
            }

            inline void ReassignAsync(std::shared_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check = false,
                                      SizeType oldVID = 0, std::function<void()> p_callback=nullptr)
            {
                auto* curJob = new ReassignAsyncJob(this, std::move(vectorContain), VID, newHeads, check, oldVID, p_callback);
                m_reassignThreadPool->add(curJob);
            }

            void ProcessAsyncReassign(std::shared_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check,
                                      SizeType oldVID, std::function<void()> p_callback);

            bool AllFinished() {return m_dispatcher->allFinished();}

            void ForceCompaction() {if (m_options.m_useKV) m_extraSearcher->ForceCompaction();}

            int getSplitTimes() {return m_splitNum;}

            int getHeadMiss() {return m_headMiss.load();}

            void UpdateStop()
            {
                m_persistentBuffer->StopPB();
                m_dispatcher->stop();
            }

            template <typename ValueType>
            void PreReassign()
            {
                //Pre Split
                bool splited = true;
                while (splited)
                {
                    splited = false;
                    for (int i = 0; i < m_extraSearcher->GetIndexSize(); i++)
                    {
                        if (m_index->ContainSample(i) && m_postingSizes[i].load() > (m_extraSearcher->GetPostingSizeLimit() * 0.8) )
                        {
                            int headID = i;
                            std::unique_lock<std::shared_timed_mutex> lock(m_rwLocks[headID]);
                            m_splitTaskNum++;
                            std::string postingList;
                            m_extraSearcher->SearchIndex(headID, postingList);

                            // reinterpret postingList to vectors and IDs
                            auto* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                            size_t vectorInfoSize = m_options.m_dim * sizeof(ValueType) + sizeof(int);
                            size_t postVectorNum = postingList.size() / vectorInfoSize;
                            COMMON::Dataset<ValueType> smallSample;  // smallSample[i] -> VID
                            std::shared_ptr<uint8_t> vectorBuffer(new uint8_t[m_options.m_dim * sizeof(ValueType) * postVectorNum], std::default_delete<uint8_t[]>());
                            std::vector<int> localIndicesInsert(postVectorNum);  // smallSample[i] = j <-> localindices[j] = i
                            std::vector<int> localIndices(postVectorNum);
                            auto vectorBuf = vectorBuffer.get();
                            size_t realVectorNum = postVectorNum;
                            int index = 0;
                            //LOG(Helper::LogLevel::LL_Info, "Scanning\n");
                            for (int j = 0; j < postVectorNum; j++)
                            {
                                uint8_t* vectorId = postingP + j * vectorInfoSize;
                                //LOG(Helper::LogLevel::LL_Info, "vector index/total:id: %d/%d:%d\n", j, m_postingSizes[headID].load(), *(reinterpret_cast<int*>(vectorId)));
                                if (CheckIdDeleted(*(reinterpret_cast<int*>(vectorId)))) {
                                    realVectorNum--;
                                } else {
                                    localIndicesInsert[index] = *(reinterpret_cast<int*>(vectorId));
                                    localIndices[index] = index;
                                    index++;
                                    memcpy(vectorBuf, vectorId + sizeof(int), m_options.m_dim * sizeof(ValueType));
                                    vectorBuf += m_options.m_dim * sizeof(ValueType);
                                }
                            }
                            // double gcEndTime = sw.getElapsedMs();
                            // m_splitGcCost += gcEndTime;
                            if (realVectorNum < (m_extraSearcher->GetPostingSizeLimit() * 0.8) )
                            {
                                postingList.clear();
                                for (int j = 0; j < realVectorNum; j++)
                                {
                                    postingList += Helper::Convert::Serialize<int>(&localIndicesInsert[j], 1);
                                    postingList += Helper::Convert::Serialize<ValueType>(vectorBuffer.get() + j * m_options.m_dim * sizeof(ValueType), m_options.m_dim);
                                }
                                m_postingSizes[headID].store(realVectorNum);
                                m_extraSearcher->AddIndex(headID, postingList);
                                // m_splitWriteBackCost += sw.getElapsedMs() - gcEndTime;
                                continue;
                            }
                            //LOG(Helper::LogLevel::LL_Info, "Resize\n");
                            localIndicesInsert.resize(realVectorNum);
                            localIndices.resize(realVectorNum);
                            smallSample.Initialize(realVectorNum, m_options.m_dim, m_index->m_iDataBlockSize, m_index->m_iDataCapacity, reinterpret_cast<ValueType*>(vectorBuffer.get()), false);

                            //LOG(Helper::LogLevel::LL_Info, "Headid: %d Sample Vector Num: %d, Real Vector Num: %d\n", headID, smallSample.R(), realVectorNum);

                            // k = 2, maybe we can change the split number, now it is fixed
                            SPTAG::COMMON::KmeansArgs<ValueType> args(2, smallSample.C(), (SizeType)localIndicesInsert.size(), 1, m_index->GetDistCalcMethod());
                            std::shuffle(localIndices.begin(), localIndices.end(), std::mt19937(std::random_device()()));
                            int numClusters = SPTAG::COMMON::KmeansClustering(smallSample, localIndices, 0, (SizeType)localIndices.size(), args);
                            if (numClusters <= 1)
                            {
                                postingList.clear();
                                float r = 0.f;
                                for (int j = 0; j < realVectorNum; j++)
                                {
                                    postingList += Helper::Convert::Serialize<int>(&localIndicesInsert[j], 1);
                                    postingList += Helper::Convert::Serialize<ValueType>(vectorBuffer.get() + j * m_options.m_dim * sizeof(ValueType), m_options.m_dim);
                                    auto dist = m_index->ComputeDistance(vectorBuffer.get() + j * m_options.m_dim * sizeof(ValueType), m_index->GetSample(headID));
                                    r = std::max<float>(r, dist);
                                }
                                m_postingSizes[headID].store(realVectorNum);
                                m_extraSearcher->AddIndex(headID, postingList);
                                continue;
                            }
                            // double clusterEndTime = sw.getElapsedMs();
                            // m_splitClusteringCost += clusterEndTime - gcEndTime;
                            splited = true;

                            long long newHeadVID = -1;
                            int first = 0;
                            std::vector<SizeType> newHeadsID(2);
                            std::vector<std::string> newPostingLists;
                            for (int k = 0; k < 2; k++) {
                                int begin, end = 0;
                                std::string postingList;
                                if (args.counts[k] == 0)	continue;

                                // LOG(Helper::LogLevel::LL_Info, "Insert new head vector\n");
                                // Notice: newHeadVID maybe an existing head vector

                                m_index->AddIndexId(smallSample[args.clusterIdx[k]], 1, m_options.m_dim, begin, end);
                                newHeadVID = begin;
                                newHeadsID.push_back(begin);
                                // LOG(Helper::LogLevel::LL_Info, "Head id: %d split into : %d, length: %d\n", headID, newHeadVID, args.counts[k]);

                                for (int j = 0; j < args.counts[k]; j++)
                                {
                                    postingList += Helper::Convert::Serialize<SizeType>(&localIndicesInsert[localIndices[first + j]], 1);
                                    postingList += Helper::Convert::Serialize<ValueType>(smallSample[localIndices[first + j]], m_options.m_dim);
                                }
                                m_extraSearcher->AddIndex(newHeadVID, postingList);
                                newPostingLists.push_back(postingList);
                                first += args.counts[k];

                                // TODO: move postingSizes to extraSearcher
                                m_postingSizes[newHeadVID] = args.counts[k];
                                m_index->AddIndexIdx(begin, end);
                            }

                            m_index->DeleteIndex(headID);
                            m_extraSearcher->DeleteIndex(headID);

                            m_postingSizes[headID] = 0;
                            lock.unlock();
                            ++m_splitNum;
                            // m_splitUpdateIndexCost += sw.getElapsedMs() - clusterEndTime;
                            //QuantifySplit(headID, newPostingLists, newHeadsID, split_order);

                            if (!m_options.m_disableReassign) ReAssign(headID, newPostingLists, newHeadsID);

                            //LOG(Helper::LogLevel::LL_Info, "After ReAssign\n");

                            //QuantifySplit(headID, newPostingLists, newHeadsID, split_order);
                        }

                        while(!AllFinished())
                        {
                            std::this_thread::sleep_for(std::chrono::milliseconds(50));
                        }
                    }
                }
            }
        };
    } // namespace SPANN
} // namespace SPTAG

#endif // _SPTAG_SPANN_INDEX_H_
