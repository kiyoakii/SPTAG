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
                std::shared_ptr<Index> m_index;
                SizeType headID;
                int appendNum;
                std::unique_ptr<std::string> appendPosting;
                std::function<void()> m_callback;
            public:
                AppendAsyncJob(std::shared_ptr<Index> m_index, SizeType headID, int appendNum, std::unique_ptr<std::string> appendPosting, std::function<void()> p_callback)
                        : m_index(m_index), headID(headID), appendNum(appendNum), appendPosting(std::move(appendPosting)), m_callback(p_callback) {}

                ~AppendAsyncJob() {}

                inline void exec() {
                    m_index->Append(headID, appendNum, appendPosting.get());
                    if (m_callback != nullptr) {
                        m_callback();
                    }
                }
            };

            class ReassignAsyncJob : public SPTAG::Helper::ThreadPool::Job
            {
            private:
                std::shared_ptr<Index> m_index;
                std::unique_ptr<std::string> vectorContain;
                SizeType VID;
                std::vector<SizeType>& newHeads;
                bool check;
                SizeType oldVID;
                std::function<void()> m_callback;
            public:
                ReassignAsyncJob(std::shared_ptr<Index> m_index,
                                 std::unique_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check,
                                 SizeType oldVID, std::function<void()> p_callback)
                        : m_index(m_index),
                          vectorContain(std::move(vectorContain)), VID(VID), newHeads(newHeads), check(check), oldVID(oldVID), m_callback(p_callback) {}

                ~ReassignAsyncJob() {}

                void exec() {
                    m_index->ProcessAsyncReassign(std::move(vectorContain), VID, newHeads, check, oldVID, std::move(m_callback));
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
            };

            class Dispatcher
            {
            private:
                std::thread t;

                std::size_t batch;
                std::atomic_bool running{false};
                std::atomic_uint32_t sentAssignment{0};

                std::shared_ptr<Index> m_index;
                std::shared_ptr<PersistentBuffer> m_persistentBuffer;
                std::shared_ptr<ThreadPool> appendThreadPool;
                std::shared_ptr<ThreadPool> reassignThreadPool;

            public:
                Dispatcher(std::shared_ptr<PersistentBuffer> pb, std::size_t batch, std::shared_ptr<ThreadPool> append, std::shared_ptr<ThreadPool> reassign)
                        : m_persistentBuffer(pb), batch(batch), appendThreadPool(append), reassignThreadPool(reassign) {}

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
            std::unique_ptr<std::shared_mutex[]> m_rwLocks;
            std::unique_ptr<std::atomic_uint32_t[]> m_postingSizes;
            std::atomic_uint64_t m_vectorNum{0};

            std::shared_ptr<IExtraSearcher> m_extraSearcher;
            std::unique_ptr<COMMON::WorkSpacePool<ExtraWorkSpace>> m_workSpacePool;

            Options m_options;

            float(*m_fComputeDistance)(const T* pX, const T* pY, DimensionType length);
            int m_iBaseSquare;

            std::shared_ptr<Dispatcher> m_dispatcher;
            std::shared_ptr<PersistentBuffer> m_persistentBuffer;
            std::unique_ptr<Helper::ThreadPool> m_threadPool;
            std::unique_ptr<ThreadPool> m_appendThreadPool;
            std::unique_ptr<ThreadPool> m_reassignThreadPool;

            COMMON::Labelset m_deletedID;
            COMMON::Labelset m_reassignedID;
            std::atomic_uint32_t m_headMiss{0};
            uint32_t m_appendTaskNum{0};
            uint32_t m_splitTaskNum{0};
            uint32_t m_split_num{0};
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

            inline SizeType GetNumSamples() const { return m_options.m_vectorSize; }
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

            ErrorCode Append(SizeType headID, int appendNum, std::string* appendPosting, SizeType oldVID);
            ErrorCode Split(const SizeType headID, int appendNum, std::string& appendPosting);
            ErrorCode ReAssign(SizeType headID, std::vector<std::string>& postingLists, std::vector<SizeType>& newHeadsID);
            void ReAssignVectors(std::map<SizeType, T*>& reAssignVectors, std::vector<SizeType>& newHeadsID, bool check=false);
            void ReAssignUpdate(std::string*, SizeType VID, std::vector<SizeType>&, bool check=false, SizeType oldVID=0);
            std::shared_ptr<COMMON::Labelset> DeletedID() {return std::make_shared<COMMON::Labelset>(m_deletedID);}

        public:
            inline void AppendAsync(SizeType headID, int appendNum, std::unique_ptr<std::string> appendPosting, std::function<void()> p_callback=nullptr)
            {
                AppendAsyncJob* curJob = new AppendAsyncJob(this, headID, appendNum, std::move(appendPosting), p_callback);
                m_appendThreadPool->add(curJob);
            }

            inline void ReassignAsync(std::unique_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check = false,
                                      SizeType oldVID = 0, std::function<void()> p_callback=nullptr)
            {
                ReassignAsyncJob* curJob = new ReassignAsyncJob(this, std::move(vectorContain), VID, newHeads, check, oldVID, p_callback);
                m_reassignThreadPool->add(curJob);
            }

            void ProcessAsyncReassign(std::unique_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check,
                                      SizeType oldVID, std::function<void()> p_callback);
        };
    } // namespace SPANN
} // namespace SPTAG

#endif // _SPTAG_SPANN_INDEX_H_
