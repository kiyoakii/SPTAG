#pragma once
#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Helper/Concurrent.h"
#include "inc/SSDServing/VectorSearch/SearchStats.h"
#include <ppl.h>

namespace SPTAG {
	namespace SSDServing {
		namespace VectorSearch {
            template <typename T>
            float CalcRecall(vector<COMMON::QueryResultSet<T>>& results, const vector<set<int>>& truth, int K)
            {
                float recall = 0;
                for (int i = 0; i < results.size(); i++)
                {
                    float thisrecall = 0;
                    for (int id : truth[i])
                    {
                        for (int j = 0; j < K; j++)
                        {
                            if (results[i].GetResult(j)->VID == id)
                            {
                                thisrecall += 1;
                                break;
                            }
                        }
                    }
                    recall += thisrecall / K;
                }
                recall /= results.size();
                return recall;
            }

            void LoadTruth(std::string truthPath, vector<set<int>>& truth, int NumQuerys, int K)
            {
                int get;
                ifstream fp(truthPath);
                if (!fp.is_open())
                {
                    fprintf(stderr, "Failed open truth file: %s\n", truthPath.c_str());
                    exit(1);
                }

                string line;
                truth.resize(NumQuerys);
                for (int i = 0; i < NumQuerys; ++i)
                {
                    truth[i].clear();
                    for (int j = 0; j < K; ++j)
                    {
                        fp >> get;
                        truth[i].insert(get);
                    }

                    std::getline(fp, line);
                }
                fp.close();
            }

            template <typename ValueType>
            void OutputResult(const std::string& p_output, vector<COMMON::QueryResultSet<ValueType>>& p_results, int p_resultNum)
            {
                if (!p_output.empty())
                {
                    std::ofstream resultOutput(p_output, std::ios::binary);
                    int32_t i32Val = static_cast<int32_t>(p_results.size());
                    resultOutput.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    i32Val = p_resultNum;
                    resultOutput.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    float fVal = 0;
                    for (size_t i = 0; i < p_results.size(); ++i)
                    {
                        for (int j = 0; j < p_resultNum; ++j)
                        {
                            i32Val = p_results[i].GetResult(j)->VID;
                            resultOutput.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                            fVal = p_results[i].GetResult(j)->Dist;
                            resultOutput.write(reinterpret_cast<char*>(&fVal), sizeof(fVal));
                        }
                    }
                }
            }

            template<typename T, typename V>
            void PrintPercentiles(FILE* p_output, const std::vector<V>& p_values, std::function<T(const V&)> p_get, const char* p_format)
            {
                double sum = 0;
                std::vector<T> collects;
                collects.reserve(p_values.size());
                for (const auto& v : p_values)
                {
                    T tmp = p_get(v);
                    sum += tmp;
                    collects.push_back(tmp);
                }

                std::sort(collects.begin(), collects.end());

                fprintf(p_output, "Avg\t50tiles\t90tiles\t95tiles\t99tiles\t99.9tiles\tMax\n");

                std::string formatStr("%.3lf");
                for (int i = 1; i < 7; ++i)
                {
                    formatStr += '\t';
                    formatStr += p_format;
                }

                formatStr += '\n';

                fprintf(p_output,
                    formatStr.c_str(),
                    sum / collects.size(),
                    collects[static_cast<size_t>(collects.size() * 0.50)],
                    collects[static_cast<size_t>(collects.size() * 0.90)],
                    collects[static_cast<size_t>(collects.size() * 0.95)],
                    collects[static_cast<size_t>(collects.size() * 0.99)],
                    collects[static_cast<size_t>(collects.size() * 0.999)],
                    collects[static_cast<size_t>(collects.size() - 1)]);
            }


            template <typename ValueType>
            void SearchSequential(SearchProcessor<ValueType>& p_searcher,
                int p_numThreads,
                vector<COMMON::QueryResultSet<ValueType>>& p_results,
                vector<SearchStats>& p_stats,
                int p_maxQueryCount,
                FILE* p_logOut)
            {
                uint32_t numQueries = std::min<uint32_t>(static_cast<uint32_t>(p_results.size()), p_maxQueryCount);
                // uint32_t subSize = (numQueries - 1) / p_numThreads + 1;

                LARGE_INTEGER frequency;
                ::QueryPerformanceFrequency(&frequency);

                LARGE_INTEGER startTime;
                ::QueryPerformanceCounter(&startTime);

                std::atomic_uint32_t processed = 0;

                fprintf(stdout, "Searching: numThread: %d, numQueries: %d.\n", p_numThreads, numQueries);
                concurrency::parallel_for(0, p_numThreads, [&](int tid)
                    {
                        // LARGE_INTEGER timePoint;
                        for (uint32_t next = processed.fetch_add(1); next < numQueries; next = processed.fetch_add(1))
                        {
                            if ((next & ((1 << 14) - 1)) == 0)
                            {
                                fprintf(stderr, "Processed %.2lf%%...\r", next * 100.0 / numQueries);
                            }

                            p_searcher.Search(p_results[next], p_stats[next], tid);
                        }
                    });

                LARGE_INTEGER finishSearch;

                ::QueryPerformanceCounter(&finishSearch);

                double sendingCost = (finishSearch.QuadPart - startTime.QuadPart) * 1.0 / frequency.QuadPart;

                fprintf(p_logOut,
                    "Finish sending in %.3lf seconds, actuallQPS is %.2lf, query count %u.\n",
                    sendingCost,
                    numQueries / sendingCost,
                    static_cast<uint32_t>(numQueries));
            }


            template <typename ValueType>
            void SearchAsync(SearchProcessor<ValueType>& p_searcher,
                uint32_t p_qps,
                vector<COMMON::QueryResultSet<ValueType>>& p_results,
                vector<SearchStats>& p_stats,
                int p_maxQueryCount,
                FILE* p_logOut)
            {
                size_t numQueries = std::min<size_t>(p_results.size(), p_maxQueryCount);

                fprintf(p_logOut, "Using Async sending with QPS setting %u\n", p_qps);

                Helper::Concurrent::WaitSignal waitFinish;
                waitFinish.Reset(static_cast<uint32_t>(numQueries));

                auto callback = [&waitFinish]()
                {
                    waitFinish.FinishOne();
                };

                std::atomic_size_t queriesSent = 0;

                LARGE_INTEGER frequency;
                ::QueryPerformanceFrequency(&frequency);

                LARGE_INTEGER startTime;
                ::QueryPerformanceCounter(&startTime);

                auto func = [&]()
                {
                    size_t index = 0;

                    while (true)
                    {
                        LARGE_INTEGER currentTime;
                        ::QueryPerformanceCounter(&currentTime);

                        float timeElapsedSec = (float)(currentTime.QuadPart - startTime.QuadPart) / frequency.QuadPart;

                        size_t targetQueries = std::min<size_t>(static_cast<size_t>(p_qps * timeElapsedSec), numQueries);

                        while (targetQueries > index)
                        {
                            index = queriesSent.fetch_add(1);

                            if (index < numQueries)
                            {
                                p_searcher.SearchAsync(p_results[index], p_stats[index], callback);

                                if ((index & ((1 << 14) - 1)) == 0)
                                {
                                    fprintf(stderr, "Sent %.2lf%%...\r", index * 100.0 / numQueries);
                                }
                            }
                            else
                            {
                                return;
                            }
                        }

                        Sleep(0);
                    }
                };

                std::thread thread1(func);
                std::thread thread2(func);
                std::thread thread3(func);

                LARGE_INTEGER finishSending;
                LARGE_INTEGER finishSearch;

                thread1.join();
                thread2.join();
                thread3.join();

                ::QueryPerformanceCounter(&finishSending);

                double sendingCost = (finishSending.QuadPart - startTime.QuadPart) * 1.0 / frequency.QuadPart;

                fprintf(p_logOut,
                    "Finish sending in %.3lf seconds, QPS setting is %u, actuallQPS is %.2lf, query count %u.\n",
                    sendingCost,
                    p_qps,
                    numQueries / sendingCost,
                    static_cast<uint32_t>(numQueries));

                waitFinish.Wait();
                ::QueryPerformanceCounter(&finishSearch);

                fprintf(p_logOut,
                    "Finish searching in %.3lf seconds.\n",
                    (finishSearch.QuadPart - startTime.QuadPart) * 1.0 / frequency.QuadPart);
            }

            template <typename ValueType>
            void Search(Options& p_opts, shared_ptr<VectorIndex> headIndex)
            {
                string queryFile = p_opts.m_queryFile;
                string outputFile = p_opts.m_searchResult;
                string truthFile = p_opts.m_truthFile;
                string warmupFile = p_opts.m_warmupFile;

                FILE* logOut = stdout;
                if (!p_opts.m_logFile.empty())
                {
                    logOut = fopen(p_opts.m_logFile.c_str(), "w");
                }

                int numThreads = p_opts.m_iNumberOfThreads;
                int asyncCallQPS = p_opts.m_qpsLimit;

                int internalResultNum = std::max<int>(p_opts.m_internalResultNum, p_opts.m_resultNum);
                int K = std::min<int>(p_opts.m_resultNum, internalResultNum);

                SearchDefault<ValueType> searcher(headIndex);
                fprintf(stderr, "Start setup index...\n");
                searcher.Setup(p_opts);

                fprintf(stderr, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, internalResultNum, asyncCallQPS > 0, p_opts);

                if (!warmupFile.empty())
                {
                    fprintf(stderr, "Start loading warmup query set...\n");
                    BasicVectorSet warmupQuerySet(warmupFile.c_str(), GetEnumValueType<ValueType>());

                    int warmupNumQueries = warmupQuerySet.Count();

                    vector<COMMON::QueryResultSet<ValueType>> warmupResults(warmupNumQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                    vector<SearchStats> warmpUpStats(warmupNumQueries);
                    for (int i = 0; i < warmupNumQueries; ++i)
                    {
                        warmupResults[i].SetTarget(reinterpret_cast<ValueType*>(warmupQuerySet.GetVector(i)));
                        warmupResults[i].Reset();
                    }

                    fprintf(stderr, "Start warmup...\n");
                    if (asyncCallQPS == 0)
                    {
                        SearchSequential(searcher, numThreads, warmupResults, warmpUpStats, p_opts.m_queryCountLimit, logOut);
                    }
                    else
                    {
                        SearchAsync(searcher, asyncCallQPS, warmupResults, warmpUpStats, p_opts.m_queryCountLimit, logOut);
                    }

                    fprintf(stderr, "\nFinish warmup...\n");
                }

                fprintf(stderr, "Start loading QuerySet...\n");
                BasicVectorSet querySet(queryFile.c_str(), GetEnumValueType<ValueType>());

                int numQueries = querySet.Count();

                vector<set<int>> truth;
                if (!truthFile.empty())
                {

                    fprintf(stderr, "Start loading TruthFile...\n");
                    LoadTruth(truthFile, truth, numQueries, K);
                }

                vector<COMMON::QueryResultSet<ValueType>> results(numQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                vector<SearchStats> stats(numQueries);
                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet.GetVector(i)));
                    results[i].Reset();
                }


                fprintf(stderr, "Start ANN Search...\n");

                LARGE_INTEGER startTime;
                LARGE_INTEGER frequency;
                ::QueryPerformanceCounter(&startTime);
                ::QueryPerformanceFrequency(&frequency);

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit, logOut);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit, logOut);
                }

                fprintf(stderr, "\nFinish ANN Search...\n");
                /*
                for (int i = 0; i < p_opts.m_queryCountLimit && i < stats.size(); ++i)
                {
                    fprintf(stderr,
                            "InQueue: %.3lf,\t OutQueue %.3lf,\t Run %.3lf,\t Exit %.3lf,\t TID %d\n",
                            (stats[i].m_searchRequestTime - startTime.QuadPart) * 1000.0 / frequency.QuadPart,
                            (stats[i].m_outQueueTS - startTime.QuadPart) * 1000.0 / frequency.QuadPart,
                            stats[i].m_totalSearchLatency,
                            (stats[i].m_exitThreadTS - startTime.QuadPart) * 1000.0 / frequency.QuadPart,
                            stats[i].m_threadID);
                }
                */

                float recall = 0;

                if (!truthFile.empty())
                {
                    recall = CalcRecall(results, truth, K);
                    fprintf(stderr, "Recall: %f\n", recall);
                }

                long long exCheckSum = 0;
                int exCheckMax = 0;
                long long exListSum = 0;
                std::for_each(stats.begin(), stats.end(), [&](const SearchStats& p_stat)
                    {
                        exCheckSum += p_stat.m_exCheck;
                        exCheckMax = std::max<int>(exCheckMax, p_stat.m_exCheck);
                        exListSum += p_stat.m_totalListElementsCount;
                    });

                fprintf(stderr,
                    "Max Ex Dist Check: %d, Average Ex Dist Check: %.2lf, Average Ex Elements Count: %.2lf.\n",
                    exCheckMax,
                    static_cast<double>(exCheckSum) / numQueries,
                    static_cast<double>(exListSum) / numQueries);

                fprintf(logOut, "\nSleep Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_sleepLatency;
                    },
                    "%.3lf");


                fprintf(logOut, "\nIn Queue Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_queueLatency;
                    },
                    "%.3lf");

                fprintf(logOut, "\nEx Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_exLatency;
                    },
                    "%.3lf");

                fprintf(logOut, "\nTotal Search Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_totalSearchLatency;
                    },
                    "%.3lf");

                fprintf(logOut, "\nTotal Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_totalLatency;
                    },
                    "%.3lf");

                fprintf(logOut, "\nTotal Disk Acess Distribution:\n");
                PrintPercentiles<int, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> int
                    {
                        return ss.m_diskAccessCount;
                    },
                    "%4d");

                fprintf(logOut, "\nTotal Async Latency 0 Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency0;
                    },
                    "%.3lf");

                fprintf(logOut, "\nTotal Async Latency 1 Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency1;
                    },
                    "%.3lf");

                fprintf(logOut, "\nTotal Async Latency 2 Distribution:\n");
                PrintPercentiles<double, SearchStats>(logOut,
                    stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency2;
                    },
                    "%.3lf");

                fprintf(logOut, "\n");

                if (!outputFile.empty())
                {
                    fprintf(stderr, "Start output to %s\n", outputFile.c_str());
                    OutputResult(outputFile, results, K);
                }

                fprintf(stdout,
                    "Recall: %f, MaxExCheck: %d, AverageExCheck: %.2lf, AverageExElements: %.2lf\n",
                    recall,
                    exCheckMax,
                    static_cast<double>(exCheckSum) / numQueries,
                    static_cast<double>(exListSum) / numQueries);

                fprintf(stdout, "\n");


                if (!p_opts.m_logFile.empty())
                {
                    fclose(logOut);
                }
            }
		}
	}
}