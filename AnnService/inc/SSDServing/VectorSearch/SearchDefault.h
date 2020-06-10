#pragma once

#ifdef _MSC_VER
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcher.h"
#else // non windows
#include "inc/SSDServing/VectorSearch/IExtraSearcherLinux.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcherLinux.h"
#endif

#include "inc/Helper/ThreadPool.h"
#include "inc/SSDServing/VectorSearch/SearchProcessor.h"
#include "inc/Core/VectorIndex.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

#include <boost/lockfree/stack.hpp>

#include <atomic>

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
                    m_tids = 0;
                }

                virtual ~SearchDefault()
                {
                }


                virtual void Setup(Options& p_config)
                {   
                    std::string vectorTranslateMap = p_config.m_vectorIDTranslate;
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

                        fprintf(stderr, "Loaded %lu Vector IDs\n", input.gcount() / sizeof(long long));
                    }

                    if (!extraFullGraphFile.empty())
                    {
                        fprintf(stderr, "Using FullGraph without cache.\n");
                        m_extraSearcher.reset(new ExtraFullGraphSearcher<ValueType>(extraFullGraphFile));
                        m_extraSearcher->Setup(p_config);
                    }
                }


                virtual void Search(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats)
                {
                    TimeUtils::StopW sw;
                    double StartingTime, EndingTime, ExEndingTime;

                    StartingTime = sw.getElapsedMs();

                    m_index->SearchIndex(p_queryResults);

                    EndingTime = sw.getElapsedMs();

                    ExtraWorkSpace* auto_ws = nullptr;
                    if (nullptr != m_extraSearcher)
                    {
                        auto_ws = GetWs();
                        auto_ws->m_postingIDs.clear();

                        for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
                        {
                            auto res = p_queryResults.GetResult(i);
                            if (res->VID != -1)
                            {
                                auto_ws->m_postingIDs.emplace_back(res->VID);
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

                        m_extraSearcher->Search(auto_ws, p_queryResults, m_index, p_stats);
                        RetWs(auto_ws);
                    }

                    ExEndingTime = sw.getElapsedMs();

                    p_stats.m_exLatency = ExEndingTime - EndingTime;
                    p_stats.m_totalSearchLatency = ExEndingTime - StartingTime;
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
                    p_stats.m_searchRequestTime = std::chrono::steady_clock::now();
                    
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
                }

                virtual shared_ptr<VectorIndex> HeadIndex() {
                    return m_index;
                }

                virtual ExtraWorkSpace* GetWs() {
                    ExtraWorkSpace* ws = nullptr;
                    if (!m_workspaces.pop(ws)) {
                        ws = new ExtraWorkSpace();
                    }

                    return ws;
                }

                virtual void RetWs(ExtraWorkSpace* ws) {
                    if (ws != nullptr)
                    {
                        m_workspaces.push(ws);
                    }
                }

            private:
                void ProcessAsyncSearch(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback)
                {
                    static thread_local int tid = m_tids.fetch_add(1);

                    std::chrono::steady_clock::time_point startPoint = std::chrono::steady_clock::now();
                    p_stats.m_queueLatency = TimeUtils::getMsInterval(p_stats.m_searchRequestTime, startPoint);

                    p_stats.m_threadID = tid;

                    static thread_local std::chrono::steady_clock::time_point m_lastQuit = startPoint;
                    p_stats.m_sleepLatency = TimeUtils::getMsInterval(m_lastQuit, startPoint);

                    Search(p_queryResults, p_stats);

                    p_stats.m_totalLatency = TimeUtils::getMsInterval(p_stats.m_searchRequestTime, std::chrono::steady_clock::now());

                    p_callback();

                    m_lastQuit = std::chrono::steady_clock::now();
                }

                shared_ptr<VectorIndex> m_index;

                std::unique_ptr<long long[]> m_vectorTranslateMap;

                std::unique_ptr<IExtraSearcher<ValueType>> m_extraSearcher;
                
                std::unique_ptr<Helper::ThreadPool> m_threadPool;

                std::atomic<std::int32_t> m_tids;

                boost::lockfree::stack<ExtraWorkSpace*> m_workspaces;
            };
        }
    }
}