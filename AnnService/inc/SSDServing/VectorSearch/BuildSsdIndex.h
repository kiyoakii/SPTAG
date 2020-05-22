#pragma once
#include <unordered_set>
#include <string>
#include <fstream>
#include <memory>
#include <vector>
#include <set>
#include <float.h>

#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/SSDServing/VectorSearch/SearchDefault.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {
            namespace Local
            {
                const std::uint16_t c_pageSize = 4096;

                struct Edge
                {
                    Edge() : headID(INT_MAX), fullID(INT_MAX), distance(FLT_MAX)
                    {
                    }

                    int headID;
                    int fullID;
                    float distance;
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
                        std::ifstream input(p_filename, std::ios::binary);
                        if (!input.is_open())
                        {
                            fprintf(stderr, "failed open VectorIDTranslate: %s\n", p_filename.c_str());
                            exit(1);
                        }

                        long long vid;
                        while (input.read(reinterpret_cast<char*>(&vid), sizeof(vid)))
                        {
                            p_set.insert(static_cast<int>(vid));
                        }

                        input.close();
                        fprintf(stderr, "Loaded %u Vector IDs\n", static_cast<uint32_t>(p_set.size()));
                    }
                    else
                    {
                        fprintf(stderr, "Not found VectorIDTranslate!\n");
                        exit(1);
                    }
                }

                void SelectPostingOffset(size_t p_spacePerVector,
                    const std::vector<std::atomic_int>& p_postingListSizes,
                    std::unique_ptr<int[]>& p_postPageNum,
                    std::unique_ptr<std::uint16_t[]>& p_postPageOffset,
                    std::vector<int>& p_postingOrderInIndex)
                {
                    p_postPageNum.reset(new int[p_postingListSizes.size()]);
                    p_postPageOffset.reset(new std::uint16_t[p_postingListSizes.size()]);

                    struct PageModWithID
                    {
                        int id;

                        std::uint16_t rest;
                    };

                    struct PageModeWithIDCmp
                    {
                        bool operator()(const PageModWithID& a, const PageModWithID& b) const
                        {
                            return a.rest == b.rest ? a.id < b.id : a.rest > b.rest;
                        }
                    };

                    std::set<PageModWithID, PageModeWithIDCmp> listRestSize;

                    p_postingOrderInIndex.clear();
                    p_postingOrderInIndex.reserve(p_postingListSizes.size());

                    PageModWithID listInfo;
                    for (size_t i = 0; i < p_postingListSizes.size(); ++i)
                    {
                        if (p_postingListSizes[i] == 0)
                        {
                            continue;
                        }

                        listInfo.id = static_cast<int>(i);
                        listInfo.rest = static_cast<std::uint16_t>((p_spacePerVector * p_postingListSizes[i]) % c_pageSize);

                        listRestSize.insert(listInfo);
                    }

                    listInfo.id = -1;

                    int currPageNum = 0;
                    std::uint16_t currOffset = 0;

                    while (!listRestSize.empty())
                    {
                        listInfo.rest = c_pageSize - currOffset;
                        auto iter = listRestSize.lower_bound(listInfo);
                        if (iter == listRestSize.end())
                        {
                            ++currPageNum;
                            currOffset = 0;
                        }
                        else
                        {
                            p_postPageNum[iter->id] = currPageNum;
                            p_postPageOffset[iter->id] = currOffset;

                            p_postingOrderInIndex.push_back(iter->id);

                            currOffset += iter->rest;
                            if (currOffset > c_pageSize)
                            {
                                fprintf(stderr, "Crossing extra pages\n");
                                exit(1);
                            }

                            if (currOffset == c_pageSize)
                            {
                                ++currPageNum;
                                currOffset = 0;
                            }

                            currPageNum += static_cast<int>((p_spacePerVector * p_postingListSizes[iter->id]) / c_pageSize);

                            listRestSize.erase(iter);
                        }
                    }

                    fprintf(stderr, "TotalPageNumbers: %d, IndexSize: %lu\n", currPageNum, static_cast<uint64_t>(currPageNum)* c_pageSize + currOffset);
                }


                void OutputSSDIndexFile(const std::string& p_outputFile,
                    size_t p_spacePerVector,
                    const std::vector<std::atomic_int>& p_postingListSizes,
                    const std::vector<Edge>& p_postingSelections,
                    const std::unique_ptr<int[]>& p_postPageNum,
                    const std::unique_ptr<std::uint16_t[]>& p_postPageOffset,
                    const std::vector<int>& p_postingOrderInIndex,
                    const BasicVectorSet& p_fullVectors)
                {
                    fprintf(stdout, "Start output...\n");

                    std::ofstream output(p_outputFile, std::ios::binary);
                    if (!output.is_open())
                    {
                        fprintf(stderr, "Failed open file %s\n", p_outputFile.c_str());
                        exit(1);
                    }

                    std::uint64_t listOffset = sizeof(int) * 4;
                    listOffset += (sizeof(int) + sizeof(std::uint16_t) + sizeof(int) + sizeof(std::uint16_t)) * p_postingListSizes.size();

                    std::unique_ptr<char[]> paddingVals(new char[c_pageSize]);
                    memset(paddingVals.get(), 0, sizeof(char) * c_pageSize);

                    std::uint64_t paddingSize = c_pageSize - (listOffset % c_pageSize);
                    if (paddingSize == c_pageSize)
                    {
                        paddingSize = 0;
                    }
                    else
                    {
                        listOffset += paddingSize;
                    }

                    // Number of lists.
                    int i32Val = static_cast<int>(p_postingListSizes.size());
                    output.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    // Number of all documents.
                    i32Val = static_cast<int>(p_fullVectors.Count());
                    output.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    // Bytes of each vector.
                    i32Val = static_cast<int>(p_fullVectors.Dimension());
                    output.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    // Page offset of list content section.
                    i32Val = static_cast<int>(listOffset / c_pageSize);
                    output.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));

                    for (int i = 0; i < p_postingListSizes.size(); ++i)
                    {
                        int pageNum = 0;
                        std::uint16_t pageOffset = 0;
                        int listEleCount = 0;
                        std::uint16_t listPageCount = 0;

                        if (p_postingListSizes[i] > 0)
                        {
                            pageNum = p_postPageNum[i];
                            pageOffset = static_cast<std::uint16_t>(p_postPageOffset[i]);
                            listEleCount = static_cast<int>(p_postingListSizes[i]);
                            listPageCount = static_cast<std::uint16_t>((p_spacePerVector * p_postingListSizes[i]) / c_pageSize);
                            if (0 != ((p_spacePerVector * p_postingListSizes[i]) % c_pageSize))
                            {
                                ++listPageCount;
                            }
                        }

                        output.write(reinterpret_cast<char*>(&pageNum), sizeof(pageNum));
                        output.write(reinterpret_cast<char*>(&pageOffset), sizeof(pageOffset));
                        output.write(reinterpret_cast<char*>(&listEleCount), sizeof(listEleCount));
                        output.write(reinterpret_cast<char*>(&listPageCount), sizeof(listPageCount));
                    }

                    if (paddingSize > 0)
                    {
                        output.write(reinterpret_cast<char*>(paddingVals.get()), paddingSize);
                    }

                    if (static_cast<uint64_t>(output.tellp()) != listOffset)
                    {
                        fprintf(stdout, "List offset not match!\n");
                        exit(1);
                    }

                    fprintf(stdout, "SubIndex Size: %lu bytes, %lu MBytes\n", listOffset, listOffset >> 20);

                    listOffset = 0;

                    std::uint64_t paddedSize = 0;
                    for (auto id : p_postingOrderInIndex)
                    {
                        std::uint64_t targetOffset = static_cast<uint64_t>(p_postPageNum[id])* c_pageSize + p_postPageOffset[id];
                        if (targetOffset < listOffset)
                        {
                            fprintf(stdout, "List offset not match, targetOffset < listOffset!\n");
                            exit(1);
                        }

                        if (targetOffset > listOffset)
                        {
                            if (targetOffset - listOffset > c_pageSize)
                            {
                                fprintf(stdout, "Padding size greater than page size!\n");
                                exit(1);
                            }

                            output.write(reinterpret_cast<char*>(paddingVals.get()), targetOffset - listOffset);

                            paddedSize += targetOffset - listOffset;

                            listOffset = targetOffset;
                        }


                        std::size_t selectIdx = std::lower_bound(p_postingSelections.begin(), p_postingSelections.end(), id, g_edgeComparer) - p_postingSelections.begin();
                        for (int j = 0; j < p_postingListSizes[id]; ++j)
                        {
                            if (p_postingSelections[selectIdx].headID != id)
                            {
                                fprintf(stdout, "Selection ID NOT MATCH\n");
                                exit(1);
                            }

                            i32Val = p_postingSelections[selectIdx++].fullID;
                            output.write(reinterpret_cast<char*>(&i32Val), sizeof(i32Val));
                            output.write(reinterpret_cast<char*>(p_fullVectors.GetVector(i32Val)), p_fullVectors.PerVectorDataSize());

                            listOffset += p_spacePerVector;
                        }
                    }

                    paddingSize = c_pageSize - (listOffset % c_pageSize);
                    if (paddingSize == c_pageSize)
                    {
                        paddingSize = 0;
                    }
                    else
                    {
                        listOffset += paddingSize;
                        paddedSize += paddingSize;
                    }

                    if (paddingSize > 0)
                    {
                        output.write(reinterpret_cast<char*>(paddingVals.get()), paddingSize);
                    }

                    output.close();

                    fprintf(stdout, "Padded Size: %lu, final total size: %lu.\n", paddedSize, listOffset);

                    fprintf(stdout, "Output done...\n");
                }
            }

            template<typename ValueType>
            void BuildSsdIndex(Options& p_opts, shared_ptr<VectorIndex> headIndex)
            {
                using namespace Local;

                TimeUtils::StopW sw;

                std::string queryFile = p_opts.m_queryFile;
                std::string outputFile = p_opts.m_ssdIndex;

                if (outputFile.empty())
                {
                    fprintf(stderr, "Output file can't be empty!\n");
                    exit(1);
                }

                int numThreads = p_opts.m_iNumberOfThreads;
                int candidateNum = p_opts.m_internalResultNum;

                std::unordered_set<int> headVectorIDS;
                LoadHeadVectorIDSet(p_opts.m_vectorIDTranslate, headVectorIDS);
                p_opts.m_vectorIDTranslate = "";

                SearchDefault<ValueType> searcher(headIndex);
                fprintf(stderr, "Start setup index...\n");
                searcher.Setup(p_opts);

                fprintf(stderr, "Setup index finish, start setup hint...\n");
                searcher.SetHint(numThreads, candidateNum, false, p_opts);

                BasicVectorSet fullVectors(queryFile.c_str(), headIndex->GetVectorValueType(), p_opts.m_iQueryDimension, p_opts.m_iQueryNumber, p_opts.m_queryFileType);

                fprintf(stderr, "Full vector loaded.\n");

                std::vector<Edge> selections;
                selections.resize(static_cast<size_t>(fullVectors.Count())* p_opts.m_replicaCount);

                std::vector<int> replicaCount(fullVectors.Count(), 0);
                std::vector<std::atomic_int> postingListSize(headIndex->GetNumSamples());
                for (auto& pls : postingListSize) pls = 0;

                fprintf(stderr, "Preparation done, start candidate searching.\n");

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
                                if (fullID >= fullVectors.Count())
                                {
                                    break;
                                }

                                if (headVectorIDS.count(fullID) > 0)
                                {
                                    continue;
                                }

                                ValueType* buffer = reinterpret_cast<ValueType*>(fullVectors.GetVector(fullID));
                                resultSet.SetTarget(buffer);
                                resultSet.Reset();

                                searcher.Search(resultSet, searchStats);

                                size_t selectionOffset = static_cast<size_t>(fullID)* p_opts.m_replicaCount;

                                BasicResult* queryResults = resultSet.GetResults();
                                for (int i = 0; i < candidateNum && replicaCount[fullID] < p_opts.m_replicaCount; ++i)
                                {
                                    if (queryResults[i].VID == -1)
                                    {
                                        continue;
                                    }

                                    // RNG Check.
                                    bool rngAccpeted = true;
                                    for (int j = 0; j < replicaCount[fullID]; ++j)
                                    {
                                        // VQANNSearch::QueryResultSet<ValueType> resultSet(NULL, candidateNum);

                                        float nnDist = headIndex->ComputeDistance(
                                            headIndex->GetSample(queryResults[i].VID),
                                            headIndex->GetSample(selections[selectionOffset + j].headID));

                                        // fprintf(stdout, "NNDist: %f Original: %f\n", nnDist, queryResults[i].Score);
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

                fprintf(stderr, "Searching replicas ended. RNG failed count: %lu\n", static_cast<uint64_t>(rngFailedCountTotal.load()));

                std::sort(selections.begin(), selections.end(), g_edgeComparer);

                int postingSizeLimit = INT_MAX;
                if (p_opts.m_postingPageLimit > 0)
                {
                    postingSizeLimit = static_cast<int>(p_opts.m_postingPageLimit * c_pageSize / (fullVectors.PerVectorDataSize() + sizeof(int)));
                }

                fprintf(stderr, "Posting size limit: %d\n", postingSizeLimit);

                {
                    std::vector<int> replicaCountDist(p_opts.m_replicaCount + 1, 0);
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        if (headVectorIDS.count(i) > 0)
                        {
                            continue;
                        }

                        ++replicaCountDist[replicaCount[i]];
                    }

                    fprintf(stderr, "Before Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        fprintf(stderr, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
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
                    std::ofstream emptyReplicaIDS("EmptyReplicaID.bin", std::ios::binary);
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        if (headVectorIDS.count(i) > 0)
                        {
                            continue;
                        }

                        ++replicaCountDist[replicaCount[i]];

                        if (replicaCount[i] < 2)
                        {
                            long long vid = i;
                            emptyReplicaIDS.write(reinterpret_cast<char*>(&vid), sizeof(vid));
                        }
                    }

                    fprintf(stderr, "After Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        fprintf(stderr, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
                    }
                }

                // VectorSize + VectorIDSize
                size_t vectorInfoSize = sizeof(ValueType) * fullVectors.Dimension() + sizeof(int);

                std::unique_ptr<int[]> postPageNum;
                std::unique_ptr<std::uint16_t[]> postPageOffset;
                std::vector<int> postingOrderInIndex;
                SelectPostingOffset(vectorInfoSize, postingListSize, postPageNum, postPageOffset, postingOrderInIndex);

                OutputSSDIndexFile(outputFile,
                    vectorInfoSize,
                    postingListSize,
                    selections,
                    postPageNum,
                    postPageOffset,
                    postingOrderInIndex,
                    fullVectors);

                double elapsedMinutes = sw.getElapsedMin();
                fprintf(stderr, "Total used time: %.2lf minutes (about %.2lf hours).\n", elapsedMinutes, elapsedMinutes / 60.0);
            }
        }
    }
}