#pragma once
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcher.h"
#include "inc/SSDServing/VectorSearch/ExtraGraphSearcher.h"
#include "inc/Helper/ThreadPool.h"
#include "inc/SSDServing/VectorSearch/SearchProcessor.h"
#include "inc/Core/VectorIndex.h"

using namespace std;

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {

            template <typename ValueType>
            class SearchDefault : public SearchProcessor<ValueType>
            {
            public:
                SearchDefault(shared_ptr<VectorIndex> p_index): m_index(p_index)
                {
                    ::QueryPerformanceFrequency(&Frequency);

                    ::QueryPerformanceCounter(&StartTime);

                    m_tids = 0;
                }

                virtual ~SearchDefault()
                {
                }


                virtual void Setup(Options& p_config)
                {   
                    std::string vectorTranslateMap = p_config.m_vectorIDTranslate;
                    std::string extraGraphFile = p_config.m_extraGraphFile;
                    std::string extraGraphVectorsetFile = p_config.m_extraGraphVectorSetFile;
                    std::string extraFullGraphFile = p_config.m_extraFullGraphFile;

                    if (!vectorTranslateMap.empty())
                    {
                        m_vectorTranslateMap.reset(new long long[m_index->GetNumSamples()]);

                        std::ifstream input(vectorTranslateMap, std::ios::binary);
                        if (!input.is_open())
                        {
                            fprintf(stderr, "failed open %s\n", vectorTranslateMap.c_str());
                            exit(1);
                        }

                        input.read(reinterpret_cast<char*>(m_vectorTranslateMap.get()), sizeof(long long)* m_index->GetNumSamples());
                        input.close();

                        fprintf(stderr, "Loaded %llu Vector IDs\n", input.gcount() / sizeof(long long));
                    }

                    if (!extraGraphFile.empty() && !extraGraphVectorsetFile.empty())
                    {
                        m_extraSearcher.reset(new ExtraGraphSearcher<ValueType>(extraGraphVectorsetFile, extraGraphFile));
                        m_extraSearcher->Setup(p_config);
                    }
                    else if (!extraFullGraphFile.empty())
                    {
                        fprintf(stderr, "Using FullGraph without cache.\n");
                        m_extraSearcher.reset(new ExtraFullGraphSearcher<ValueType>(extraFullGraphFile));
                        m_extraSearcher->Setup(p_config);
                    }
                }


                virtual void Search(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, int p_threadID)
                {
                    LARGE_INTEGER StartingTime, EndingTime, ExEndingTime;

                    QueryPerformanceCounter(&StartingTime);

                    m_index->SearchIndex(p_queryResults);

                    QueryPerformanceCounter(&EndingTime);

                    if (nullptr != m_extraSearcher)
                    {
                        m_extraWorkSpaces[p_threadID]->m_postingIDs.clear();

                        for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
                        {
                            auto res = p_queryResults.GetResult(i);
                            if (res->VID != -1)
                            {
                                m_extraWorkSpaces[p_threadID]->m_postingIDs.emplace_back(res->VID);
                            }
                        }
                    }

                    if (m_vectorTranslateMap != nullptr)
                    {
                        for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
                        {
                            auto res = p_queryResults.GetResult(i);
                            if (res->VID != -1)
                            {
                                res->VID = static_cast<int>(m_vectorTranslateMap[res->VID]);
                            }
                        }
                    }

                    if (nullptr != m_extraSearcher)
                    {
                        p_queryResults.Reverse();

                        m_extraSearcher->Search(m_extraWorkSpaces[p_threadID], p_queryResults, m_index, p_stats);
                    }

                    QueryPerformanceCounter(&ExEndingTime);

                    p_stats.m_exLatency = (ExEndingTime.QuadPart - EndingTime.QuadPart) * 1000.0 / Frequency.QuadPart;
                    p_stats.m_totalSearchLatency = (ExEndingTime.QuadPart - StartingTime.QuadPart) * 1000.0 / Frequency.QuadPart;
                    p_stats.m_totalLatency = p_stats.m_totalSearchLatency;
                }

                class SearchAsyncJob : public SPTAG::Helper::ThreadPool::Job
                {
                private:
                    SearchDefault* m_processor;
                    COMMON::QueryResultSet<ValueType>& m_queryResults;
                    SearchStats& m_stats;
                    std::function<void()> m_callback;
                public:
                    SearchAsyncJob(SearchDefault* p_processor,
                        COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback)
                        : m_processor(p_processor),
                        m_queryResults(p_queryResults), m_stats(p_stats), m_callback(p_callback) {}

                    ~SearchAsyncJob() {}

                    void exec() {
                        m_processor->ProcessAsyncSearch(m_queryResults, m_stats, std::move(m_callback));
                    }
                };

                virtual void SearchAsync(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback)
                {
                    LARGE_INTEGER t;
                    QueryPerformanceCounter(&t);
                    p_stats.m_searchRequestTime = t.QuadPart;
                    
                    SearchAsyncJob* curJob = new SearchAsyncJob(this, p_queryResults, p_stats, p_callback);

                    m_threadPool->add(curJob);
                }


                virtual void SetHint(int p_threadNum, int p_resultNum, bool p_asyncCall, const Options& p_opts)
                {
                    fprintf(stderr, "ThreadNum: %d, ResultNum: %d, AsyncCall: %d\n", p_threadNum, p_resultNum, p_asyncCall ? 1 : 0);

                    if (p_asyncCall)
                    {
                        m_threadPool.reset(new Helper::ThreadPool());
                        m_threadPool->init(p_threadNum);
                    }

                    if (nullptr != m_extraSearcher)
                    {
                        m_extraWorkSpaces.resize(p_threadNum, nullptr);
                        for (size_t i = 0; i < p_threadNum; i++)
                        {
                            m_extraWorkSpaces[i] = new ExtraWorkSpace();
                        }
                    }
                }

                virtual shared_ptr<VectorIndex> HeadIndex() {
                    return m_index;
                }

            private:
                void ProcessAsyncSearch(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback)
                {
                    static thread_local int tid = m_tids.fetch_add(1);

                    LARGE_INTEGER t;
                    QueryPerformanceCounter(&t);
                    p_stats.m_queueLatency = (t.QuadPart - p_stats.m_searchRequestTime) * 1000.0 / Frequency.QuadPart;

                    p_stats.m_threadID = tid;
                    p_stats.m_outQueueTS = t.QuadPart;

                    static thread_local LARGE_INTEGER m_lastQuit{ t };
                    p_stats.m_sleepLatency = (t.QuadPart - m_lastQuit.QuadPart) * 1000.0 / Frequency.QuadPart;

                    Search(p_queryResults, p_stats, tid);

                    QueryPerformanceCounter(&t);

                    p_stats.m_totalLatency = (t.QuadPart - p_stats.m_searchRequestTime) * 1000.0 / Frequency.QuadPart;

                    p_callback();

                    QueryPerformanceCounter(&m_lastQuit);

                    p_stats.m_exitThreadTS = m_lastQuit.QuadPart;
                }

                shared_ptr<VectorIndex> m_index;

                std::unique_ptr<long long[]> m_vectorTranslateMap;

                std::unique_ptr<IExtraSearcher<ValueType>> m_extraSearcher;
                
                std::vector<ExtraWorkSpace*> m_extraWorkSpaces;

                LARGE_INTEGER Frequency;

                LARGE_INTEGER StartTime;

                std::unique_ptr<Helper::ThreadPool> m_threadPool;

                std::atomic_int32_t m_tids;
            };
        }
    }
}