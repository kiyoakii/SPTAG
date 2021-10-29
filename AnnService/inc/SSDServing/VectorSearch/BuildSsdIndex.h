#pragma once
#include <unordered_set>
#include <string>
#include <memory>
#include <vector>
#include <set>
#include <float.h>

#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"
#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/SSDServing/VectorSearch/SearchDefault.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Helper/VectorSetReader.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {
            namespace Local
            {
                const std::uint16_t c_pageSize = 4096;

                struct Edge
                {
                    Edge() : headID(INT_MAX), fullID(INT_MAX), distance(FLT_MAX), order(0)
                    {
                    }

                    int headID;
                    int fullID;
                    float distance;
					char order;
                };

                struct EdgeCompare
                {
                    bool operator()(const Edge& a, int b) const
                    {
                        return a.headID < b;
                    };

                    bool operator()(int a, const Edge& b) const
                    {
                        return a < b.headID;
                    };

                    bool operator()(const Edge& a, const Edge& b) const
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
                } g_edgeComparer;

                void LoadHeadVectorIDSet(const std::string& p_filename, std::unordered_set<int>& p_set)
                {
                    if (!p_filename.empty())
                    {
                        auto ptr = SPTAG::f_createIO();
                        if (ptr == nullptr || !ptr->Initialize(p_filename.c_str(), std::ios::binary | std::ios::in)) {
                            LOG(Helper::LogLevel::LL_Error, "failed open VectorIDTranslate: %s\n", p_filename.c_str());
                            exit(1);
                        }

                        long long vid;
                        while (ptr->ReadBinary(sizeof(vid), reinterpret_cast<char*>(&vid)) == sizeof(vid))
                        {
                            p_set.insert(static_cast<int>(vid));
                        }
                        LOG(Helper::LogLevel::LL_Info, "Loaded %u Vector IDs\n", static_cast<uint32_t>(p_set.size()));
                    }
                    else
                    {
                        LOG(Helper::LogLevel::LL_Error, "Not found VectorIDTranslate!\n");
                        exit(1);
                    }
                }
            }

            template<typename ValueType>
            void BuildSsdIndex(Options& p_opts)
            {
                using namespace Local;

                TimeUtils::StopW sw;

                std::string outputFile = COMMON_OPTS.m_ssdIndex;

                if (outputFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Error, "Output file can't be empty!\n");
                    exit(1);
                }

                int numThreads = p_opts.m_iNumberOfThreads;
                int candidateNum = p_opts.m_internalResultNum;

                std::unordered_set<int> headVectorIDS;
                LoadHeadVectorIDSet(COMMON_OPTS.m_headIDFile, headVectorIDS);

                SearchDefault<ValueType> searcher;
                LOG(Helper::LogLevel::LL_Info, "Start setup index...\n");
                searcher.Setup(p_opts);

                LOG(Helper::LogLevel::LL_Info, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, candidateNum, false, p_opts);

                std::shared_ptr<Helper::ReaderOptions> vectorOptions(new Helper::ReaderOptions(COMMON_OPTS.m_valueType, COMMON_OPTS.m_dim, COMMON_OPTS.m_vectorType, COMMON_OPTS.m_vectorDelimiter));
                auto vectorReader = Helper::VectorSetReader::CreateInstance(vectorOptions);
                if (ErrorCode::Success != vectorReader->LoadFile(COMMON_OPTS.m_vectorPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read vector file.\n");
                    exit(1);
                }
                auto fullVectors = vectorReader->GetVectorSet();
                if (COMMON_OPTS.m_distCalcMethod == DistCalcMethod::Cosine) fullVectors->Normalize(p_opts.m_iNumberOfThreads);

                LOG(Helper::LogLevel::LL_Info, "Full vector loaded.\n");

                std::vector<Edge> selections(static_cast<size_t>(fullVectors->Count())* p_opts.m_replicaCount);

                std::vector<int> replicaCount(fullVectors->Count(), 0);
                std::vector<std::atomic_int> postingListSize(searcher.HeadIndex()->GetNumSamples());
                for (auto& pls : postingListSize) pls = 0;

                LOG(Helper::LogLevel::LL_Info, "Preparation done, start candidate searching.\n");

                std::vector<std::thread> threads;
                threads.reserve(numThreads);

                std::atomic_int nextFullID(0);
                std::atomic_size_t rngFailedCountTotal(0);

                for (int tid = 0; tid < numThreads; ++tid)
                {
                    threads.emplace_back([&, tid]()
                        {
                            COMMON::QueryResultSet<ValueType> resultSet(NULL, candidateNum);
                            SearchStats searchStats;

                            size_t rngFailedCount = 0;

                            while (true)
                            {
                                int fullID = nextFullID.fetch_add(1);
                                if (fullID >= fullVectors->Count())
                                {
                                    break;
                                }

                                if (!COMMON_OPTS.m_addHeadToPost && headVectorIDS.count(fullID) > 0)
                                {
                                    continue;
                                }

                                ValueType* buffer = reinterpret_cast<ValueType*>(fullVectors->GetVector(fullID));
                                resultSet.SetTarget(buffer);
                                resultSet.Reset();

                                searcher.Search(resultSet, searchStats);

                                size_t selectionOffset = static_cast<size_t>(fullID)* p_opts.m_replicaCount;

                                BasicResult* queryResults = resultSet.GetResults();
                                for (int i = 0; i < candidateNum && replicaCount[fullID] < p_opts.m_replicaCount; ++i)
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

                                        float nnDist = searcher.HeadIndex()->ComputeDistance(
                                            searcher.HeadIndex()->GetSample(queryResults[i].VID),
                                            searcher.HeadIndex()->GetSample(selections[selectionOffset + j].headID));

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

                for (int tid = 0; tid < numThreads; ++tid)
                {
                    threads[tid].join();
                }

                LOG(Helper::LogLevel::LL_Info, "Searching replicas ended. RNG failed count: %llu\n", static_cast<uint64_t>(rngFailedCountTotal.load()));

                std::sort(selections.begin(), selections.end(), g_edgeComparer);

                int postingSizeLimit = INT_MAX;
                if (p_opts.m_postingPageLimit > 0)
                {
                    postingSizeLimit = static_cast<int>(p_opts.m_postingPageLimit * c_pageSize / (fullVectors->PerVectorDataSize() + sizeof(int)));
                }
                if (COMMON_OPTS.m_addHeadToPost)
                {
                    postingSizeLimit += 1;
                }

                LOG(Helper::LogLevel::LL_Info, "Posting size limit: %d\n", postingSizeLimit);

                {
                    std::vector<int> replicaCountDist(p_opts.m_replicaCount + 1, 0);
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        if (!COMMON_OPTS.m_addHeadToPost && headVectorIDS.count(i) > 0)
                        {
                            continue;
                        }

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
                    if (postingListSize[i] <= postingSizeLimit)
                    {
                        continue;
                    }

                    std::size_t selectIdx = std::lower_bound(selections.begin(), selections.end(), i, g_edgeComparer) - selections.begin();
                    for (size_t dropID = postingSizeLimit; dropID < postingListSize[i]; ++dropID)
                    {
                        int fullID = selections[selectIdx + dropID].fullID;
                        --replicaCount[fullID];
                    }

                    postingListSize[i] = postingSizeLimit;
                }

                if (p_opts.m_outputEmptyReplicaID)
                {
                    std::vector<int> replicaCountDist(p_opts.m_replicaCount + 1, 0);
                    auto ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize("EmptyReplicaID.bin", std::ios::binary | std::ios::out)) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to create EmptyReplicaID.bin!\n");
                        exit(1);
                    }
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        if (!COMMON_OPTS.m_addHeadToPost && headVectorIDS.count(i) > 0)
                        {
                            continue;
                        }

                        ++replicaCountDist[replicaCount[i]];

                        if (replicaCount[i] < 2)
                        {
                            long long vid = i;
                            if (ptr->WriteBinary(sizeof(vid), reinterpret_cast<char*>(&vid)) != sizeof(vid)) {
                                LOG(Helper::LogLevel::LL_Error, "Failt to write EmptyReplicaID.bin!");
                                exit(1);
                            }
                        }
                    }

                    LOG(Helper::LogLevel::LL_Info, "After Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        LOG(Helper::LogLevel::LL_Info, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
                    }
                }

                std::string postinglist;
                for (int id = 0; id < postingListSize.size(); id++) 
                {
                    postinglist.resize(0);
                    postinglist.clear();
                    std::size_t selectIdx = std::lower_bound(selections.begin(), selections.end(), id, g_edgeComparer)
                                            - selections.begin();
                    for (int j = 0; j < postingListSize[id]; ++j) {
                        if (selections[selectIdx].headID != id) {
                            LOG(Helper::LogLevel::LL_Error, "Selection ID NOT MATCH\n");
                            exit(1);
                        }
                        int fullID = selections[selectIdx++].fullID;
                        size_t dim = fullVectors->Dimension();
                        // First Vector ID, then Vector
                        postinglist += Helper::Serialize<int>(&fullID, 1);
                        postinglist += Helper::Serialize<ValueType>(fullVectors->GetVector(fullID), dim);
                    }

                    db->Put(WriteOptions(), Helper::Serialize<int>(&id, 1), postinglist);
                }
                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(COMMON_OPTS.m_ssdIndexInfo.c_str(), std::ios::binary | std::ios::out))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", COMMON_OPTS.m_ssdIndexInfo.c_str());
                    exit(1);
                }
                //Number of all documents.
                int i32Val = static_cast<int>(fullVectors->Count());
                if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to write SSDIndexInfo File!");
                    exit(1);
                }
                //Number of postings
                i32Val = static_cast<int>(postingListSize.size());
                if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to write SSDIndexInfo File!");
                    exit(1);
                }
                for(int id = 0; id < postingListSize.size(); id++)
                {
                    i32Val = postingListSize[id];
                    if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to write SSDIndexInfo File!");
                        exit(1);
                    }
                }
                i32Val = static_cast<int>(fullVectors->Count());
                COMMON::Labelset m_deletedID;
                m_deletedID.Initialize(i32Val);
                m_deletedID.Save(COMMON_OPTS.m_deleteID);
                
                double elapsedMinutes = sw.getElapsedMin();
                LOG(Helper::LogLevel::LL_Info, "Total used time: %.2lf minutes (about %.2lf hours).\n", elapsedMinutes, elapsedMinutes / 60.0);
            }
        }
    }
}
