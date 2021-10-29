#pragma once

#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"
#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Helper/Concurrent.h"
#include "inc/SSDServing/VectorSearch/SearchStats.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"
#include "inttypes.h"

namespace SPTAG {
	namespace SSDServing {
		namespace VectorSearch {
            template <typename T>
            float CalcRecall(std::vector<COMMON::QueryResultSet<T>>& results, const std::vector<std::set<int>>& truth, int K)
            {
                float recall = 0;
                for (int i = 0; i < results.size(); i++)
                {
                    for (int id : truth[i])
                    {
                        //LOG(Helper::LogLevel::LL_Info, "truth id: %d ", id);
                        for (int j = 0; j < K; j++)
                        {
                            //LOG(Helper::LogLevel::LL_Info, "vector id %d: %d", j, results[i].GetResult(j)->VID);
                            if (results[i].GetResult(j)->VID == id)
                            {
                                recall++;
                                break;
                            }
                        }
                        //LOG(Helper::LogLevel::LL_Info, "\n");
                    }
                }
                return static_cast<float>(recall)/static_cast<float>(results.size() * K);
            }

            void LoadTruthTXT(std::string truthPath, std::vector<std::set<int>>& truth, int K, SizeType p_iTruthNumber)
            {
                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(truthPath.c_str(), std::ios::in)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed open truth file: %s\n", truthPath.c_str());
                    exit(1);
                }
                std::size_t lineBufferSize = 20;
                std::unique_ptr<char[]> currentLine(new char[lineBufferSize]);
                truth.clear();
                truth.resize(p_iTruthNumber);
                for (int i = 0; i < p_iTruthNumber; ++i)
                {
                    truth[i].clear();
                    for (int j = 0; j < K; ++j)
                    {
                        if (ptr->ReadString(lineBufferSize, currentLine, ' ') == 0) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to read truth file!\n");
                            exit(1);
                        }
                        truth[i].insert(std::atoi(currentLine.get()));
                    }
                    if (ptr->ReadString(lineBufferSize, currentLine, '\n') == 0) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to read truth file!\n");
                        exit(1);
                    }
                }
            }

            void LoadTruthXVEC(std::string truthPath, std::vector<std::set<int>>& truth, int K, SizeType p_iTruthNumber)
            {
                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(truthPath.c_str(), std::ios::in | std::ios::binary)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed open truth file: %s\n", truthPath.c_str());
                    exit(1);
                }

                DimensionType dim = K;
                std::vector<int> temp_vec(K);
                truth.clear();
                truth.resize(p_iTruthNumber);
                for (size_t i = 0; i < p_iTruthNumber; i++) {
                    if (ptr->ReadBinary(4, (char*)&dim) != 4 || dim < K) {
                        LOG(Helper::LogLevel::LL_Error, "Error: Xvec file %s has No.%" PRId64 " vector whose dims are fewer than expected. Expected: %" PRId32 ", Fact: %" PRId32 "\n", truthPath.c_str(), i, K, dim);
                        exit(1);
                    }
                    if (dim > K) temp_vec.resize(dim);
                    if (ptr->ReadBinary(dim * 4, (char*)temp_vec.data()) != dim * 4) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to read truth file!\n");
                        exit(1);
                    }
                    truth[i].insert(temp_vec.begin(), temp_vec.begin() + K);
                }
            }

            void LoadTruthDefault(std::string truthPath, std::vector<std::set<int>>& truth, int K, SizeType p_iTruthNumber) {
                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(truthPath.c_str(), std::ios::in | std::ios::binary)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed open truth file: %s\n", truthPath.c_str());
                    exit(1);
                }

                int row, column;
                if (ptr->ReadBinary(4, (char*)&row) != 4 || ptr->ReadBinary(4, (char*)&column) != 4) {
                    LOG(Helper::LogLevel::LL_Error, "Fail to read truth file!\n");
                    exit(1);
                }
                truth.clear();
                truth.resize(row);
                std::vector<int> vec(column);
                for (size_t i = 0; i < row; i++)
                {
                    if (ptr->ReadBinary(4 * column, (char*)vec.data()) != 4 * column) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to read truth file!\n");
                        exit(1);
                    }
                    //LOG(Helper::LogLevel::LL_Info, "read truth size %d: %d\n", i, vec.size());
                    //LOG(Helper::LogLevel::LL_Info, "need to insert %d truth\n", K);
                    truth[i].insert(vec.begin(), vec.begin() + K);
                    //LOG(Helper::LogLevel::LL_Info, "inserted truth size %d: %d\n", i, truth[i].size());
                }
            }

            void LoadTruth(std::string truthPath, std::vector<std::set<int>>& truth, int NumQuerys, int K)
            {
                if (COMMON_OPTS.m_truthType == TruthFileType::TXT)
                {
                    LoadTruthTXT(truthPath, truth, K, NumQuerys);
                } 
                else if (COMMON_OPTS.m_truthType == TruthFileType::XVEC)
                {
                    LoadTruthXVEC(truthPath, truth, K, NumQuerys);
                }
                else if (COMMON_OPTS.m_truthType == TruthFileType::DEFAULT) {
                    LoadTruthDefault(truthPath, truth, K, NumQuerys);
                }
                else
                {
                    LOG(Helper::LogLevel::LL_Error, "TruthFileType Unsupported.\n");
                    exit(1);
                }
            }

            template <typename ValueType>
            void OutputResult(const std::string& p_output, std::vector<COMMON::QueryResultSet<ValueType>>& p_results, int p_resultNum)
            {
                if (!p_output.empty())
                {
                    auto ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize(p_output.c_str(), std::ios::binary | std::ios::out)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed create file: %s\n", p_output.c_str());
                        exit(1);
                    }
                    int32_t i32Val = static_cast<int32_t>(p_results.size());
                    if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to write result file!\n");
                        exit(1);
                    }
                    i32Val = p_resultNum;
                    if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to write result file!\n");
                        exit(1);
                    }

                    float fVal = 0;
                    for (size_t i = 0; i < p_results.size(); ++i)
                    {
                        for (int j = 0; j < p_resultNum; ++j)
                        {
                            i32Val = p_results[i].GetResult(j)->VID;
                            if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to write result file!\n");
                                exit(1);
                            }

                            fVal = p_results[i].GetResult(j)->Dist;
                            if (ptr->WriteBinary(sizeof(fVal), reinterpret_cast<char*>(&fVal)) != sizeof(fVal)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to write result file!\n");
                                exit(1);
                            }
                        }
                    }
                }
            }

            template<typename T, typename V>
            void PrintPercentiles(const std::vector<V>& p_values, std::function<T(const V&)> p_get, const char* p_format)
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

                LOG(Helper::LogLevel::LL_Info, "Avg\t50tiles\t90tiles\t95tiles\t99tiles\t99.9tiles\tMax\n");

                std::string formatStr("%.3lf");
                for (int i = 1; i < 7; ++i)
                {
                    formatStr += '\t';
                    formatStr += p_format;
                }

                formatStr += '\n';

                LOG(Helper::LogLevel::LL_Info,
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
            void SearchSequential(SearchDefault<ValueType>& p_searcher,
                int p_numThreads,
                std::vector<COMMON::QueryResultSet<ValueType>>& p_results,
                std::vector<SearchStats>& p_stats,
                int p_maxQueryCount)
            {
                int numQueries = min(static_cast<int>(p_results.size()), p_maxQueryCount);

                TimeUtils::StopW sw;

                LOG(Helper::LogLevel::LL_Info, "Searching: numThread: %d, numQueries: %d.\n", p_numThreads, numQueries);
#pragma omp parallel for num_threads(p_numThreads)
                for (int next = 0; next < numQueries; ++next)
                {
                    if ((next & ((1 << 14) - 1)) == 0)
                    {
                        LOG(Helper::LogLevel::LL_Info, "Processed %.2lf%%...\n", next * 100.0 / numQueries);
                    }

                    p_searcher.Search(p_results[next], p_stats[next]);
                }

                double sendingCost = sw.getElapsedSec();

                LOG(Helper::LogLevel::LL_Info,
                    "Finish sending in %.3lf seconds, actuallQPS is %.2lf, query count %u.\n",
                    sendingCost,
                    numQueries / sendingCost,
                    static_cast<uint32_t>(numQueries));
            }


            template <typename ValueType>
            void SearchAsync(SearchDefault<ValueType>& p_searcher,
                uint32_t p_qps,
                std::vector<COMMON::QueryResultSet<ValueType>>& p_results,
                std::vector<SearchStats>& p_stats,
                int p_maxQueryCount)
            {
                size_t numQueries = std::min<size_t>(p_results.size(), p_maxQueryCount);

                LOG(Helper::LogLevel::LL_Info, "Using Async sending with QPS setting %u\n", p_qps);

                Helper::Concurrent::WaitSignal waitFinish;
                waitFinish.Reset(static_cast<uint32_t>(numQueries));

                auto callback = [&waitFinish]()
                {
                    waitFinish.FinishOne();
                };

                std::atomic_size_t queriesSent(0);

                TimeUtils::SteadClock::time_point startTime = TimeUtils::SteadClock::now();

                auto func = [&]()
                {
                    size_t index = 0;

                    while (true)
                    {
                        TimeUtils::SteadClock::time_point currentTime = TimeUtils::SteadClock::now();

                        double timeElapsedSec = TimeUtils::getMsInterval(startTime, currentTime);

                        size_t targetQueries = std::min<size_t>(static_cast<size_t>(p_qps * timeElapsedSec), numQueries);

                        while (targetQueries > index)
                        {
                            index = queriesSent.fetch_add(1);

                            if (index < numQueries)
                            {
                                p_searcher.SearchAsync(p_results[index], p_stats[index], callback);

                                if ((index & ((1 << 14) - 1)) == 0)
                                {
                                    LOG(Helper::LogLevel::LL_Info, "Sent %.2lf%%...\n", index * 100.0 / numQueries);
                                }
                            }
                            else
                            {
                                return;
                            }
                        }
                    }
                };

                std::thread thread1(func);
                std::thread thread2(func);
                std::thread thread3(func);

                thread1.join();
                thread2.join();
                thread3.join();

                TimeUtils::SteadClock::time_point finishSending = TimeUtils::SteadClock::now();
                double sendingCost = TimeUtils::getSecInterval(startTime, finishSending);

                LOG(Helper::LogLevel::LL_Info,
                    "Finish sending in %.3lf seconds, QPS setting is %u, actuallQPS is %.2lf, query count %u.\n",
                    sendingCost,
                    p_qps,
                    numQueries / sendingCost,
                    static_cast<uint32_t>(numQueries));

                waitFinish.Wait();

                TimeUtils::SteadClock::time_point finishSearch = TimeUtils::SteadClock::now();
                double searchCost = TimeUtils::getSecInterval(startTime, finishSearch);

                LOG(Helper::LogLevel::LL_Info,
                    "Finish searching in %.3lf seconds.\n",
                    searchCost);
            }

            template <typename ValueType>
            void TestUpdateAcc(Options& p_opts)
            {
                std::string outputFile = p_opts.m_searchResult;

                if (!p_opts.m_logFile.empty())
                {
                    SPTAG::g_pLogger.reset(new Helper::FileLogger(Helper::LogLevel::LL_Info, p_opts.m_logFile.c_str()));
                }
                int numThreads = p_opts.m_iNumberOfThreads;
                int asyncCallQPS = p_opts.m_qpsLimit;

                int internalResultNum = std::max<int>(p_opts.m_internalResultNum, 64);
                int K = std::min<int>(p_opts.m_resultNum, internalResultNum);

                SearchDefault<ValueType> searcher;
                LOG(Helper::LogLevel::LL_Info, "Start setup index...\n");
                searcher.Setup(p_opts);

                LOG(Helper::LogLevel::LL_Info, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, internalResultNum, asyncCallQPS > 0, p_opts);

                searcher.LoadDeleteID(COMMON_OPTS.m_deleteID);

                LOG(Helper::LogLevel::LL_Info, "Start loading QuerySet...\n");
                std::shared_ptr<Helper::ReaderOptions> queryOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_queryType, COMMON_OPTS.m_queryDelimiter));
                auto queryReader = Helper::VectorSetReader::CreateInstance(queryOptions);
                LOG(Helper::LogLevel::LL_Info, "Start Reading QuerySet...\n");
                if (ErrorCode::Success != queryReader->LoadFile(COMMON_OPTS.m_queryPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
                    exit(1);
                }
                LOG(Helper::LogLevel::LL_Info, "Start Getting QuerySet...\n");
                auto querySet = queryReader->GetVectorSet();
                int numQueries = querySet->Count();

                std::vector<COMMON::QueryResultSet<ValueType>> results(numQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                std::vector<SearchStats> stats(numQueries);
                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }


                LOG(Helper::LogLevel::LL_Info, "Start First ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "\nFinish First ANN Search...\n");
                
                LOG(Helper::LogLevel::LL_Info, "Start Checking Update\n");

                LOG(Helper::LogLevel::LL_Info, "Testing Delete\n");

                int deletenum = std::min<int>(K, 5);
                std::set<int>deletedPoints;

                for (int i = 0; i < results.size(); i++)
                {
                    for (int j = 0; j < deletenum; j++)
                    {
                        deletedPoints.insert(results[i].GetResult(j)->VID);
                        searcher.Delete(results[i].GetResult(j)->VID);
                    }
                }

                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }

                LOG(Helper::LogLevel::LL_Info, "Start Second ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "\nFinish Second ANN Search...\n");

                LOG(Helper::LogLevel::LL_Info, "\nTesting Delete Result...\n");

                for (int i = 0; i < results.size(); i++)
                {
                    for (int j = 0; j < results[i].GetResultNum(); j++)
                    {
                        if (deletedPoints.count(results[i].GetResult(j)->VID))
                        {
                            LOG(Helper::LogLevel::LL_Error, "Error! Find deleted point %d, query: %d\n", results[i].GetResult(j)->VID, i);
                            exit(1);
                        }
                    }
                }

                LOG(Helper::LogLevel::LL_Info, "\nTesting Insert...\n");

                int insertnum = std::min<int>(numQueries, 10000);
                for (int i = 0; i < insertnum; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }
                for (int i = 0; i < insertnum; i++)
                {
                    searcher.Insert(results[i], stats[i]);
                }

                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }

                LOG(Helper::LogLevel::LL_Info, "Start Thrid ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "Finish Thrid ANN Search...\n");
                
                LOG(Helper::LogLevel::LL_Info, "Testing Insert Result...\n");
                bool findOne = false;
                for (int i = 0; i < insertnum; i++)
                {
                    //LOG(Helper::LogLevel::LL_Info, "%d : result id: ", i);
                    for (int j = 0; j < results[i].GetResultNum(); j++)
                    {
                        //LOG(Helper::LogLevel::LL_Info, "%d ", results[i].GetResult(j)->VID);
                        if (results[i].GetResult(j)->Dist < 1e-6)
                        {
                            findOne = true;
                        }
                    }
                    //LOG(Helper::LogLevel::LL_Info, "\n");
                }
                if (findOne)
                {
                    LOG(Helper::LogLevel::LL_Info, "Insert Successfully!\n");
                }
                else
                {
                    LOG(Helper::LogLevel::LL_Error, "Insert Error!\n");
                    exit(1);
                }
            }
            
            template <typename ValueType>
            void TestUpdateSta(Options& p_opts)
            {
                LOG(Helper::LogLevel::LL_Info, "Start testing update stablity\n");

                std::string truthFile = COMMON_OPTS.m_truthPath;
                std::string outputFile = p_opts.m_searchResult;

                if (!p_opts.m_logFile.empty())
                {
                    SPTAG::g_pLogger.reset(new Helper::FileLogger(Helper::LogLevel::LL_Info, p_opts.m_logFile.c_str()));
                }
                int numThreads = p_opts.m_iNumberOfThreads;
                int asyncCallQPS = p_opts.m_qpsLimit;

                int internalResultNum = std::max<int>(p_opts.m_internalResultNum, 64);
                int K = std::min<int>(p_opts.m_resultNum, internalResultNum);
                int cycle = 1;

                SearchDefault<ValueType> searcher;
                LOG(Helper::LogLevel::LL_Info, "Start setup index...\n");
                searcher.Setup(p_opts);

                LOG(Helper::LogLevel::LL_Info, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, internalResultNum, asyncCallQPS > 0, p_opts);

                searcher.LoadDeleteID(COMMON_OPTS.m_deleteID);

                LOG(Helper::LogLevel::LL_Info, "Start loading VectorSet...\n");
                std::shared_ptr<Helper::ReaderOptions> vectorOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_vectorType, COMMON_OPTS.m_vectorDelimiter));
                auto vectorReader = Helper::VectorSetReader::CreateInstance(vectorOptions);
                if (ErrorCode::Success != vectorReader->LoadFile(COMMON_OPTS.m_vectorPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read vector file.\n");
                    exit(1);
                }
                auto fullVectors = vectorReader->GetVectorSet();

                int updateVectorNum = fullVectors->Count() * p_opts.m_indexSize;

                LOG(Helper::LogLevel::LL_Info, "Start loading QuerySet...\n");
                std::shared_ptr<Helper::ReaderOptions> queryOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_queryType, COMMON_OPTS.m_queryDelimiter));
                auto queryReader = Helper::VectorSetReader::CreateInstance(queryOptions);
                if (ErrorCode::Success != queryReader->LoadFile(COMMON_OPTS.m_queryPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
                    exit(1);
                }

                std::vector<std::set<int>> truth;

                auto querySet = queryReader->GetVectorSet();
                int numQueries = querySet->Count();

                if (!truthFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Info, "Start loading TruthFile...\n");
                    LoadTruth(truthFile, truth, numQueries, K);
                }

                std::vector<COMMON::QueryResultSet<ValueType>> results(numQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                std::vector<SearchStats> statsPrev(numQueries);
                std::vector<SearchStats> stats(numQueries);
                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }


                LOG(Helper::LogLevel::LL_Info, "Start First ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, statsPrev, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, statsPrev, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "Finish First ANN Search...\n");

                float recall = 0;
                std::vector<int> recallPrev;
                recallPrev.clear();
                recallPrev.resize(numQueries);
                for (int i = 0; i < numQueries; i++)
                {
                    recallPrev[i] = 0;
                }

                std::vector<int> recallCur;
                recallCur.clear();
                recallCur.resize(numQueries);
                for (int i = 0; i < numQueries; i++)
                {
                    recallCur[i] = 0;
                }

                if (!truthFile.empty())
                {
                    recall = 0;
                    for (int k = 0; k < results.size(); k++)
                    {
                        for (int id : truth[k])
                        {
                            for (int j = 0; j < K; j++)
                            {
                                if (results[k].GetResult(j)->VID == id)
                                {
                                    recallPrev[k]++;
                                    recall++;
                                    break;
                                }
                            }
                        }
                    }
                    recall = static_cast<float>(recall)/static_cast<float>(results.size() * K);
                    LOG(Helper::LogLevel::LL_Info, "Recall: %f\n", recall);
                }
                //std::string storeFilePrev = p_opts.m_SSDVectorDistPath;
                //storeFilePrev += "_prev";
                //searcher.CalDBDist(storeFilePrev);
                /*
                if (!outputFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Info, "Start output to %s\n", outputFile.c_str());
                    OutputResult(outputFile, results, K);
                }
                */

                LOG(Helper::LogLevel::LL_Info, "Index delete/re-insert index size : %lf\n", p_opts.m_indexSize);

                std::vector<COMMON::QueryResultSet<ValueType>> insertResult(updateVectorNum, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                
                std::vector<SizeType> updateIndice;
                std::vector<SizeType> indices;

                indices.clear();
                indices.resize(fullVectors->Count());
                for (int i = 0; i < fullVectors->Count(); i++) 
                {
                    indices[i] = i;
                }


                int vectorNum = fullVectors->Count();
                int totalNum = fullVectors->Count();
                int section = totalNum / numThreads;
                int selectionNum = updateVectorNum / numThreads;
                for (int i = 0; i < cycle; i++)
                {

                    for (int k = 0; k < numQueries; k++)
                    {
                        recallCur[k] = 0;
                    }

                    if (p_opts.m_randomDisabled)
                    {
                        auto ptr = SPTAG::f_createIO();
                        if (ptr == nullptr || !ptr->Initialize(p_opts.m_insertVectorsPath.c_str(), std::ios::binary | std::ios::in))
                        {
                            LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", p_opts.m_insertVectorsPath.c_str());
                            exit(1);
                        }
                        updateIndice.resize(updateVectorNum);
                        if (ptr->ReadBinary(sizeof(int) * updateVectorNum, (char*)updateIndice.data()) != sizeof(int) * updateVectorNum) {
                            LOG(Helper::LogLevel::LL_Error, "Failed to read update vector number!\n");
                            exit(1);
                        }
                    }
                    else
                    {
                        updateIndice.clear();
                        updateIndice.resize(updateVectorNum);
                        LOG(Helper::LogLevel::LL_Info, "cycle %d: selecting update vector\n", i);
                        omp_set_num_threads(numThreads);
                        #pragma omp parallel
                            {
                                int tid = omp_get_thread_num();
                                int lowerBound = tid * section;
                                int upperBound = (tid == numThreads-1) ? totalNum : (tid + 1) * section;
                                int selectBegin = tid * selectionNum;
                                int select = (tid == numThreads-1) ? (updateVectorNum- selectionNum * (numThreads - 1)) : selectionNum;
                                LOG(Helper::LogLevel::LL_Info, "tid: %d, select range: (%d, %d), select indice range: (%d,%d), select num: %d\n", tid, lowerBound, upperBound, selectBegin, selectBegin+select, select);
                                std::vector<int> tempIndice;
                                tempIndice.resize(select);
                                for (int k = 0; k < select; k++)
                                {
                                    SizeType randid = COMMON::Utils::rand(upperBound, lowerBound);
                                    while (std::count(tempIndice.begin(), tempIndice.end(), randid)) {
                                        randid = COMMON::Utils::rand(upperBound, lowerBound);
                                    }
                                    updateIndice[k + selectBegin] = randid;
                                    tempIndice[k] = randid;
                                }
                            }
                    }
                    /*
                    for (int k = 0; k < updateVectorNum; k++) 
                    {
                        SizeType randid = COMMON::Utils::rand(fullVectors->Count(), 0);
                        while (std::count(updateIndice.begin(), updateIndice.end(), randid)) {
                            randid = COMMON::Utils::rand(fullVectors->Count(), 0);
                        }
                        updateIndice[k] = randid;
                    }
                    */
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: update vector\n", i);
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: delete vector\n", i);
                    for (int k = 0; k < updateVectorNum; k++) 
                    {
                        searcher.Delete(indices[updateIndice[k]]);
                    }
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: prepare vector\n", i);
                    for (int k = 0; k < updateVectorNum; k++)
                    {
                        //LOG(Helper::LogLevel::LL_Info, "cycle %d: insert %d vector %d \n", i, k, updateIndice[k]);
                        insertResult[k].SetTarget(reinterpret_cast<ValueType*>(fullVectors->GetVector(updateIndice[k])));
                        insertResult[k].Reset();
                    }
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: insert vector\n", i);
                    for (int k = 0; k < updateVectorNum; k++) 
                    {
                        if ((k+1) % 10000 == 0) LOG(Helper::LogLevel::LL_Info, "cycle %d: inserted %d vectors\n", i, k+1);
                        searcher.Insert(insertResult[k], stats[k]);
                        indices[updateIndice[k]] = vectorNum++;
                    }
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: after %d insertion, head vectors split %d times\n", i, updateVectorNum, searcher.getSplitNum());
                    searcher.setSplitZero();

                    LOG(Helper::LogLevel::LL_Info, "cycle %d: setup search\n", i);
                    for (int j = 0; j < numQueries; ++j)
                    {
                        results[j].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(j)));
                        results[j].Reset();
                    }
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: begin searching\n", i);
                    if (asyncCallQPS == 0)
                    {
                        SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                    }
                    else
                    {
                        SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                    }
                    LOG(Helper::LogLevel::LL_Info, "cycle %d: begin calclating recall\n", i);
                    

                    int numLow = 0;
                    int totalCountPrev = 0;
                    int totalCountCur = 0;
                    if (!truthFile.empty())
                    {
                        recall = 0;
                        for (int k = 0; k < results.size(); k++)
                        {
                            for (int id : truth[k])
                            {
                                for (int j = 0; j < K; j++)
                                {
                                    if (results[k].GetResult(j)->VID == indices[id])
                                    {
                                        recallCur[k]++;
                                        recall++;
                                        break;
                                    }
                                }
                            }
                            if (recallCur[k] < recallPrev[k])
                            {
                                //LOG(Helper::LogLevel::LL_Info, "cycle %d: query %d get prev result:\n", i, k);
                                numLow++;
                                int countPrev = 0;
                                for (std::map<int,int>::iterator it=statsPrev[k].m_headAndPostingSize.begin(); it!=statsPrev[k].m_headAndPostingSize.end(); ++it)
                                {
                                    countPrev += it->second;
                                    //LOG(Helper::LogLevel::LL_Info, "headid: %d, size: %d ", it->first, it->second);
                                }
                                //LOG(Helper::LogLevel::LL_Info, "\n");
                                //LOG(Helper::LogLevel::LL_Info, "cycle %d: query %d get cur result:\n", i, k);
                                int countCur = 0;
                                for (std::map<int,int>::iterator it=stats[k].m_headAndPostingSize.begin(); it!=stats[k].m_headAndPostingSize.end(); ++it)
                                {
                                    countCur += it->second;
                                    //LOG(Helper::LogLevel::LL_Info, "headid: %d, size: %d ", it->first, it->second);
                                }
                                //LOG(Helper::LogLevel::LL_Info, "cycle %d: query:%d get\n", i, k);
                                //LOG(Helper::LogLevel::LL_Info, "prev hit: %d checked size: %d\n", recallPrev[k], countPrev);
                                //LOG(Helper::LogLevel::LL_Info, "curr hit: %d checked size: %d\n", recallCur[k], countCur);
                                //LOG(Helper::LogLevel::LL_Info, "\n");
                                totalCountPrev += countPrev;
                                totalCountCur += countCur;
                                /*
                                if (recallCur[k] + 3 < recallPrev[k])
                                {
                                    LOG(Helper::LogLevel::LL_Info, "prev hit: %d checked size: %d\n", recallPrev[k], countPrev);
                                    LOG(Helper::LogLevel::LL_Info, "curr hit: %d checked size: %d\n", recallCur[k], countCur);
                                    LOG(Helper::LogLevel::LL_Info, "cycle %d: query %d get prev result:\n", i, k);
                                    for (std::map<int,int>::iterator it=statsPrev[k].m_headAndPostingSize.begin(); it!=statsPrev[k].m_headAndPostingSize.end(); ++it)
                                    {
                                        LOG(Helper::LogLevel::LL_Info, "headid: %d, size: %d ", it->first, it->second);
                                    }
                                    LOG(Helper::LogLevel::LL_Info, "\n");
                                    LOG(Helper::LogLevel::LL_Info, "cycle %d: query %d get cur result:\n", i, k);
                                    for (std::map<int,int>::iterator it=stats[k].m_headAndPostingSize.begin(); it!=stats[k].m_headAndPostingSize.end(); ++it)
                                    {
                                        LOG(Helper::LogLevel::LL_Info, "headid: %d, size: %d ", it->first, it->second);
                                    }
                                    LOG(Helper::LogLevel::LL_Info, "\n");
                                }
                                */
                            }
                        }
                        recall = static_cast<float>(recall)/static_cast<float>(results.size() * K);
                        LOG(Helper::LogLevel::LL_Info, "cycle %d: Recall: %f\n", i, recall);
                        LOG(Helper::LogLevel::LL_Info, "Checked Compare(Recall Low): Prev: %d Curr: %d\n", totalCountPrev/numLow, totalCountCur/numLow);
                        totalCountPrev = 0;
                        totalCountCur = 0;
                        for (int k = 0; k < numQueries; k++)
                        {
                            for (std::map<int,int>::iterator it=statsPrev[k].m_headAndPostingSize.begin(); it!=statsPrev[k].m_headAndPostingSize.end(); ++it)
                            {
                                totalCountPrev += it->second;
                            }
                            for (std::map<int,int>::iterator it=stats[k].m_headAndPostingSize.begin(); it!=stats[k].m_headAndPostingSize.end(); ++it)
                            {
                                totalCountCur += it->second;
                            }
                        }
                        LOG(Helper::LogLevel::LL_Info, "Checked Compare(Total): Prev: %d Curr: %d\n", totalCountPrev/numQueries, totalCountCur/numQueries);
                    }

                    //std::string storeFileCurr = p_opts.m_SSDVectorDistPath;
                    //storeFileCurr += "_curr";
                    //searcher.CalDBDist(storeFileCurr);

                    if (!p_opts.m_randomDisabled)
                    {
                        //store the inserted vector ID
                        auto ptr = SPTAG::f_createIO();
                        if (ptr == nullptr || !ptr->Initialize(p_opts.m_insertVectorsPath.c_str(), std::ios::binary | std::ios::out))
                        {
                            LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", p_opts.m_insertVectorsPath.c_str());
                            exit(1);
                        }
                        if (ptr->WriteBinary(sizeof(int) * updateVectorNum, (char*)updateIndice.data()) != sizeof(int) * updateVectorNum) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to store inserted vector ID");
                            exit(1);
                        }
                    }

                    /*
                    auto sptr = SPTAG::f_createIO();
                    if (sptr == nullptr || !sptr->Initialize(p_opts.m_headDistPostingnum.c_str(), std::ios::binary | std::ios::out))
                    {
                        LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", p_opts.m_headDistPostingnum.c_str());
                        exit(1);
                    }
                    for (int k = 0; k < numQueries; k++)
                    {
                        if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&recallPrev[k])) != sizeof(int)) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to store recallPrev, query number: %d\n", k);
                            exit(1);
                        }

                        if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&recallCur[k])) != sizeof(int)) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to store recallCurr, query number: %d\n", k);
                            exit(1);
                        }

                        int postingSize = statsPrev[k].m_headAndDist.size();
                        if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&postingSize)) != sizeof(int)) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev, query number: %d\n", k);
                            exit(1);
                        }
                        for (std::map<int,float>::iterator it=statsPrev[k].m_headAndDist.begin(); it!=statsPrev[k].m_headAndDist.end(); ++it)
                        {
                            int headID = it->first;
                            if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&headID)) != sizeof(int)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev headID, query number: %d\n", k);
                                exit(1);
                            }

                            float headDist = it->second;
                            if (sptr->WriteBinary(sizeof(float), reinterpret_cast<char*>(&headDist)) != sizeof(float)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev Dist, query number: %d\n", k);
                                exit(1);
                            }

                            int headSize = statsPrev[k].m_headAndPostingSize[it->first];
                            if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&headSize)) != sizeof(int)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev posting size, query number: %d\n", k);
                                exit(1);
                            }
                        }
                        postingSize = stats[k].m_headAndDist.size();
                        if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&postingSize)) != sizeof(int)) {
                            LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev, query number: %d\n", k);
                            exit(1);
                        }
                        for (std::map<int,float>::iterator it=stats[k].m_headAndDist.begin(); it!=stats[k].m_headAndDist.end(); ++it)
                        {
                            int headID = it->first;
                            if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&headID)) != sizeof(int)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev headID, query number: %d\n", k);
                                exit(1);
                            }

                            float headDist = it->second;
                            if (sptr->WriteBinary(sizeof(float), reinterpret_cast<char*>(&headDist)) != sizeof(float)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev Dist, query number: %d\n", k);
                                exit(1);
                            }

                            int headSize = stats[k].m_headAndPostingSize[it->first];
                            if (sptr->WriteBinary(sizeof(int), reinterpret_cast<char*>(&headSize)) != sizeof(int)) {
                                LOG(Helper::LogLevel::LL_Error, "Fail to store statsPrev posting size, query number: %d\n", k);
                                exit(1);
                            }
                        }
                    }
                    if (!outputFile.empty())
                    {
                        LOG(Helper::LogLevel::LL_Info, "Start output to %s %d\n", outputFile.c_str(), i);
                        OutputResult(outputFile + Helper::Serialize<int>(&i, 1), results, K);
                    }
                    */
                }
            }

            std::string GetTruthFileName(std::string& truthFilePrefix, int vectorCount)
            {
                std::string fileName(truthFilePrefix);
                fileName += "-";
                if (vectorCount < 1000)
                {
                    fileName += std::to_string(vectorCount);
                } 
                else if (vectorCount < 1000000)
                {
                    fileName += std::to_string(vectorCount/1000);
                    fileName += "k";
                }
                else if (vectorCount < 1000000000)
                {
                    fileName += std::to_string(vectorCount/1000000);
                    fileName += "M";
                }
                else
                {
                    fileName += std::to_string(vectorCount/1000000000);
                    fileName += "B";
                }
                return fileName;
            }

            template <typename ValueType>
            void TestIncremental(Options& p_opts)
            {
                LOG(Helper::LogLevel::LL_Info, "Start incremental test\n");

                std::string truthFilePrefix = p_opts.m_truthFilePrefix;
                int step = p_opts.m_step;
                if (step == 0)
                {
                    LOG(Helper::LogLevel::LL_Error, "Incremental Test Error, Need to set step.\n");
                    exit(1);
                }

                if (!p_opts.m_logFile.empty())
                {
                    SPTAG::g_pLogger.reset(new Helper::FileLogger(Helper::LogLevel::LL_Info, p_opts.m_logFile.c_str()));
                }
                int numThreads = p_opts.m_iNumberOfThreads;
                int asyncCallQPS = p_opts.m_qpsLimit;

                int internalResultNum = std::max<int>(p_opts.m_internalResultNum, 64);
                int K = std::min<int>(p_opts.m_resultNum, internalResultNum);

                SearchDefault<ValueType> searcher;
                LOG(Helper::LogLevel::LL_Info, "Start setup index...\n");
                searcher.Setup(p_opts);

                LOG(Helper::LogLevel::LL_Info, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, internalResultNum, asyncCallQPS > 0, p_opts);

                searcher.LoadDeleteID(COMMON_OPTS.m_deleteID);

                LOG(Helper::LogLevel::LL_Info, "Start loading VectorSet...\n");
                std::shared_ptr<Helper::ReaderOptions> vectorOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_vectorType, COMMON_OPTS.m_vectorDelimiter));
                auto vectorReader = Helper::VectorSetReader::CreateInstance(vectorOptions);
                if (ErrorCode::Success != vectorReader->LoadFile(COMMON_OPTS.m_extraVectorPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read extra vector file.\n");
                    exit(1);
                }
                auto extraVectors = vectorReader->GetVectorSet();

                LOG(Helper::LogLevel::LL_Info, "Start loading QuerySet...\n");
                std::shared_ptr<Helper::ReaderOptions> queryOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_queryType, COMMON_OPTS.m_queryDelimiter));
                auto queryReader = Helper::VectorSetReader::CreateInstance(queryOptions);
                if (ErrorCode::Success != queryReader->LoadFile(COMMON_OPTS.m_queryPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
                    exit(1);
                }

                std::vector<std::set<int>> truth;

                auto querySet = queryReader->GetVectorSet();
                int numQueries = querySet->Count();
                int insertCount = extraVectors->Count();
                int curCount = searcher.getVecNum();

                std::string truthfile;

                LOG(Helper::LogLevel::LL_Info, "Start loading TruthFile...\n");
                LoadTruth(GetTruthFileName(truthFilePrefix, curCount), truth, numQueries, K);

                std::vector<COMMON::QueryResultSet<ValueType>> results(numQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                std::vector<SearchStats> stats(numQueries);
                std::vector<COMMON::QueryResultSet<ValueType>> insertResults(insertCount, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }

                LOG(Helper::LogLevel::LL_Info, "Start ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "\nFinish ANN Search...\n");

                float recall = 0;

                recall = CalcRecall(results, truth, K);
                LOG(Helper::LogLevel::LL_Info, "Recall: %f\n", recall);

                for (int i = 0; i < insertCount; i++)
                {
                    insertResults[i].SetTarget(reinterpret_cast<ValueType*>(extraVectors->GetVector(i)));
                    insertResults[i].Reset();
                }
                for (int i = 0; i < insertCount; i++)
                {
                    searcher.Insert(insertResults[i], stats[i]);
                    if ((i+1) % 10000 == 0) LOG(Helper::LogLevel::LL_Info, "inserted %d vectors\n", i+1);
                    if ((i+1) % step == 0)
                    {
                        curCount += step;
                        LOG(Helper::LogLevel::LL_Info, "Total Vector num %d \n", curCount);
                        LOG(Helper::LogLevel::LL_Info, "Start Searching\n");
                        for (int i = 0; i < numQueries; ++i)
                        {
                            results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                            results[i].Reset();
                        }
                        if (asyncCallQPS == 0)
                        {
                            SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                        }
                        else
                        {
                            SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                        }
                        LOG(Helper::LogLevel::LL_Info, "Start loading TruthFile...\n");
                        LoadTruth(GetTruthFileName(truthFilePrefix, curCount), truth, numQueries, K);
                        recall = CalcRecall(results, truth, K);
                        LOG(Helper::LogLevel::LL_Info, "Recall: %f\n", recall);
                    }
                }

                searcher.Rebuild();

                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }

                LOG(Helper::LogLevel::LL_Info, "Start ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "\nFinish ANN Search...\n");

                LOG(Helper::LogLevel::LL_Info, "Start loading TruthFile...\n");
                LoadTruth(GetTruthFileName(truthFilePrefix, curCount), truth, numQueries, K);
                recall = CalcRecall(results, truth, K);
                LOG(Helper::LogLevel::LL_Info, "Recall: %f\n", recall);
            }

            template <typename ValueType>
            void Search(Options& p_opts)
            {
                if (COMMON_OPTS.m_testUpdateAcc)
                {
                    TestUpdateAcc<ValueType>(p_opts);
                    return;
                } else if (COMMON_OPTS.m_testUpdateSta)
                {
                    TestUpdateSta<ValueType>(p_opts);
                    return;
                } else if (COMMON_OPTS.m_testInc)
                {
                    TestIncremental<ValueType>(p_opts);
                    return;
                }
                std::string outputFile = p_opts.m_searchResult;
                std::string truthFile = COMMON_OPTS.m_truthPath;
                std::string warmupFile = COMMON_OPTS.m_warmupPath;

                if (!p_opts.m_logFile.empty())
                {
                    SPTAG::g_pLogger.reset(new Helper::FileLogger(Helper::LogLevel::LL_Info, p_opts.m_logFile.c_str()));
                }
                int numThreads = p_opts.m_iNumberOfThreads;
                int asyncCallQPS = p_opts.m_qpsLimit;

                int internalResultNum = std::max<int>(p_opts.m_internalResultNum, p_opts.m_resultNum);
                int K = std::min<int>(p_opts.m_resultNum, internalResultNum);

                SearchDefault<ValueType> searcher;
                LOG(Helper::LogLevel::LL_Info, "Start setup index...\n");
                searcher.Setup(p_opts);

                LOG(Helper::LogLevel::LL_Info, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, internalResultNum, asyncCallQPS > 0, p_opts);
                searcher.LoadDeleteID(COMMON_OPTS.m_deleteID);

                if (!warmupFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Info, "Start loading warmup query set...\n");
                    std::shared_ptr<Helper::ReaderOptions> queryOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_warmupType, COMMON_OPTS.m_warmupDelimiter));
                    auto queryReader = Helper::VectorSetReader::CreateInstance(queryOptions);
                    if (ErrorCode::Success != queryReader->LoadFile(COMMON_OPTS.m_warmupPath))
                    {
                        LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
                        exit(1);
                    }
                    auto warmupQuerySet = queryReader->GetVectorSet();
                    int warmupNumQueries = warmupQuerySet->Count();

                    std::vector<COMMON::QueryResultSet<ValueType>> warmupResults(warmupNumQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                    std::vector<SearchStats> warmpUpStats(warmupNumQueries);
                    for (int i = 0; i < warmupNumQueries; ++i)
                    {
                        warmupResults[i].SetTarget(reinterpret_cast<ValueType*>(warmupQuerySet->GetVector(i)));
                        warmupResults[i].Reset();
                    }

                    LOG(Helper::LogLevel::LL_Info, "Start warmup...\n");
                    if (asyncCallQPS == 0)
                    {
                        SearchSequential(searcher, numThreads, warmupResults, warmpUpStats, p_opts.m_queryCountLimit);
                    }
                    else
                    {
                        SearchAsync(searcher, asyncCallQPS, warmupResults, warmpUpStats, p_opts.m_queryCountLimit);
                    }

                    LOG(Helper::LogLevel::LL_Info, "\nFinish warmup...\n");
                }

                LOG(Helper::LogLevel::LL_Info, "Start loading QuerySet...\n");
                std::shared_ptr<Helper::ReaderOptions> queryOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_queryType, COMMON_OPTS.m_queryDelimiter));
                auto queryReader = Helper::VectorSetReader::CreateInstance(queryOptions);
                if (ErrorCode::Success != queryReader->LoadFile(COMMON_OPTS.m_queryPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
                    exit(1);
                }
                auto querySet = queryReader->GetVectorSet();
                int numQueries = querySet->Count();

                std::vector<std::set<int>> truth;
                if (!truthFile.empty())
                {

                    LOG(Helper::LogLevel::LL_Info, "Start loading TruthFile...\n");
                    LoadTruth(truthFile, truth, numQueries, K);
                }
                /*
                for (int i = 0; i < numQueries; i++)
                {
                    LOG(Helper::LogLevel::LL_Info, "query size %d: %d\n", i, truth[i].size());
                }
                */

                std::vector<COMMON::QueryResultSet<ValueType>> results(numQueries, COMMON::QueryResultSet<ValueType>(NULL, internalResultNum));
                std::vector<SearchStats> stats(numQueries);
                for (int i = 0; i < numQueries; ++i)
                {
                    results[i].SetTarget(reinterpret_cast<ValueType*>(querySet->GetVector(i)));
                    results[i].Reset();
                }


                LOG(Helper::LogLevel::LL_Info, "Start ANN Search...\n");

                if (asyncCallQPS == 0)
                {
                    SearchSequential(searcher, numThreads, results, stats, p_opts.m_queryCountLimit);
                }
                else
                {
                    SearchAsync(searcher, asyncCallQPS, results, stats, p_opts.m_queryCountLimit);
                }

                LOG(Helper::LogLevel::LL_Info, "\nFinish ANN Search...\n");

                float recall = 0;

                if (!truthFile.empty())
                {
                    recall = CalcRecall(results, truth, K);
                    LOG(Helper::LogLevel::LL_Info, "Recall: %f\n", recall);
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

                LOG(Helper::LogLevel::LL_Info,
                    "Max Ex Dist Check: %d, Average Ex Dist Check: %.2lf, Average Ex Elements Count: %.2lf.\n",
                    exCheckMax,
                    static_cast<double>(exCheckSum) / numQueries,
                    static_cast<double>(exListSum) / numQueries);

                LOG(Helper::LogLevel::LL_Info, "\nSleep Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_sleepLatency;
                    },
                    "%.3lf");


                LOG(Helper::LogLevel::LL_Info, "\nIn Queue Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_queueLatency;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nEx Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_exLatency;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Search Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_totalSearchLatency;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Latency Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_totalLatency;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Disk Acess Distribution:\n");
                PrintPercentiles<int, SearchStats>(stats,
                    [](const SearchStats& ss) -> int
                    {
                        return ss.m_diskAccessCount;
                    },
                    "%4d");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Async Latency 0 Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency0;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Async Latency 1 Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency1;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\nTotal Async Latency 2 Distribution:\n");
                PrintPercentiles<double, SearchStats>(stats,
                    [](const SearchStats& ss) -> double
                    {
                        return ss.m_asyncLatency2;
                    },
                    "%.3lf");

                LOG(Helper::LogLevel::LL_Info, "\n");

                if (!outputFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Info, "Start output to %s\n", outputFile.c_str());
                    OutputResult(outputFile, results, K);
                }

                LOG(Helper::LogLevel::LL_Info,
                    "Recall: %f, MaxExCheck: %d, AverageExCheck: %.2lf, AverageExElements: %.2lf\n",
                    recall,
                    exCheckMax,
                    static_cast<double>(exCheckSum) / numQueries,
                    static_cast<double>(exListSum) / numQueries);

                LOG(Helper::LogLevel::LL_Info, "\n");
            }
		}
	}
}