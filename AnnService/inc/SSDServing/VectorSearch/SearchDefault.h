#pragma once

#ifdef _MSC_VER
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcher.h"
#else // non windows
#include "inc/SSDServing/VectorSearch/IExtraSearcherLinux.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcherLinux.h"
#endif

#include "inc/Core/BKT/Index.h"
#include "inc/Core/Common/BKTreeUpdate.h"
#include "inc/Core/Common/Dataset.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Core/VectorIndex.h"
#include "inc/Helper/ThreadPool.h"
#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"
#include "inc/SSDServing/IndexBuildManager/Utils.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

#include <boost/lockfree/stack.hpp>

#include <atomic>
#include <shared_mutex>

namespace SPTAG {
	namespace SSDServing {
		namespace VectorSearch {
			// LARGE_INTEGER g_systemPerfFreq;
			const std::uint16_t p_pageSize = 4096;

			struct EdgeInsert
            {
                EdgeInsert() : headID(INT_MAX), fullID(INT_MAX), distance(FLT_MAX), order(0)
                {
                }

                int headID;
                int fullID;
                float distance;
				char order;
            };

			template <typename ValueType>
			class SearchDefault
			{
			public:
				SearchDefault()
					:m_workspaces(128)
				{
					m_tids = 0;
					m_replicaCount = 4;
					//QueryPerformanceFrequency(&g_systemPerfFreq);
				}

				~SearchDefault()
				{
					ExtraWorkSpace* context;
					while (m_workspaces.pop(context))
					{
						delete context;
					}
				}

				void LoadHeadIndex(Options& p_opts) {
					LOG(Helper::LogLevel::LL_Info, "Start loading head index. \n");

					if (VectorIndex::LoadIndex(COMMON_OPTS.m_headIndexFolder, m_index) != ErrorCode::Success) {
						LOG(Helper::LogLevel::LL_Error, "ERROR: Cannot Load index files!\n");
						exit(1);
					}

					m_index->SetParameter("NumberOfThreads", std::to_string(p_opts.m_iNumberOfThreads));
					m_index->SetParameter("MaxCheck", std::to_string(p_opts.m_maxCheck));

					if (!p_opts.m_headConfig.empty())
					{
						Helper::IniReader iniReader;
						if (iniReader.LoadIniFile(p_opts.m_headConfig) != ErrorCode::Success) {
							LOG(Helper::LogLevel::LL_Error, "ERROR of loading head index config: %s\n", p_opts.m_headConfig.c_str());
							exit(1);
						}

						for (const auto& iter : iniReader.GetParameters("Index"))
						{
							m_index->SetParameter(iter.first.c_str(), iter.second.c_str());
						}
					}

					LOG(Helper::LogLevel::LL_Info, "End loading head index. \n");
				}

				void LoadVectorIdsSSDIndex(std::string vectorTranslateMap, std::string extraFullGraphFile)
				{
					if (vectorTranslateMap.empty()) {
						LOG(Helper::LogLevel::LL_Error, "Config error: VectorTranlateMap Empty for Searching SSD vectors.\n");
						exit(1);
					}

					/*
					if (extraFullGraphFile.empty()) {
						LOG(Helper::LogLevel::LL_Error, "Config error: SsdIndex empty for Searching SSD vectors.\n");
						exit(1);
					}
					*/
					//m_vectorTranslateMap.reserve(m_index->GetNumSamples() * 2);

					std::unique_ptr<long long[]> tempMap;
					tempMap.reset(new long long[m_index->GetNumSamples()]);

					auto ptr = f_createIO();
					if (ptr == nullptr || !ptr->Initialize(vectorTranslateMap.c_str(), std::ios::binary | std::ios::in)) {
						LOG(Helper::LogLevel::LL_Error, "Failed open %s\n", vectorTranslateMap.c_str());
						exit(1);
					}
					if (ptr->ReadBinary(sizeof(long long) * m_index->GetNumSamples(), reinterpret_cast<char*>(tempMap.get())) != sizeof(long long) * m_index->GetNumSamples()) {
						LOG(Helper::LogLevel::LL_Error, "Failed to read vectorTanslateMap!\n");
						exit(1);
					}

					m_vectorTranslateMap.Initialize(m_index->GetNumSamples(), 1, tempMap.get(), false);
					m_vectorTranslateMap.SetName("Head");

					LOG(Helper::LogLevel::LL_Info, "Using FullGraph without cache.\n");
					std::ifstream input(COMMON_OPTS.m_ssdIndexInfo, std::ios::binary);
                    if (!input.is_open())
                    {
                        fprintf(stderr, "Failed to open file: %s\n", COMMON_OPTS.m_ssdIndexInfo.c_str());
                        exit(1);
                    }

                    input.read(reinterpret_cast<char*>(&m_vectornum), sizeof(m_vectornum));

					LOG(Helper::LogLevel::LL_Info, "Current vector num: %d.\n", m_vectornum);

					input.close();

					m_extraSearcher.reset(new ExtraFullGraphSearcher<ValueType>(extraFullGraphFile));
				}

				void CheckHeadIndexType() {
					SPTAG::VectorValueType v1 = m_index->GetVectorValueType(), v2 = GetEnumValueType<ValueType>();
					if (v1 != v2) {
						LOG(Helper::LogLevel::LL_Error, "Head index and vectors don't have the same value types, which are %s %s\n",
							SPTAG::Helper::Convert::ConvertToString(v1).c_str(),
							SPTAG::Helper::Convert::ConvertToString(v2).c_str()
						);
						exit(1);
					}
				}

				void LoadIndex4ANNIndexTestTool(const std::string& p_config,
					const std::vector<ByteArray>& p_indexBlobs,
					std::string& vectorTranslateMap,
					std::string& extraFullGraphFile)
				{
					if (VectorIndex::LoadIndex(p_config, p_indexBlobs, m_index) != SPTAG::ErrorCode::Success) {
						LOG(Helper::LogLevel::LL_Error, "LoadIndex error in LoadIndex4ANNIndexTestTool.\n");
						exit(1);
					}
					CheckHeadIndexType();
					LoadVectorIdsSSDIndex(vectorTranslateMap, extraFullGraphFile);
				}

				void LoadDeleteID(std::string m_deleteIDFilename)
				{
					
					m_deletedID.Load(m_deleteIDFilename);
				}

				void Setup(Options& p_config)
				{
					LoadHeadIndex(p_config);
					CheckHeadIndexType();
					if (!p_config.m_buildSsdIndex)
					{
						LoadVectorIdsSSDIndex(COMMON_OPTS.m_headIDFile, COMMON_OPTS.m_ssdIndex);
					}
					m_clearHead = !p_config.m_buildSsdIndex;
				}


				ErrorCode InsertPostingList(SPTAG::COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, SizeType VID) 
				{
					//fetch postingList & update

					//check the size of postingList
					//if there are oversize postingList
					//clustering, update headvector postingList
					//add new vecror
					if (COMMON_OPTS.m_indexAlgoType == IndexAlgoType::BKT) {
						int replicaCount = 0;
						BasicResult* queryResults = p_queryResults.GetResults();
						std::vector<EdgeInsert> selections(static_cast<size_t>(m_replicaCount));
						for (int i = 0; i < m_internalResultNum && replicaCount < m_replicaCount; ++i)
                        {
                            if (queryResults[i].VID == -1)
                            {
                                break;
                            }
							//LOG(Helper::LogLevel::LL_Info, "Head Vector: %d\n", queryResults[i].VID);

                            // RNG Check.
                            bool rngAccpeted = true;
                            for (int j = 0; j < replicaCount; ++j)
                            {
								float nnDist = m_index->ComputeDistance(
                                m_index->GetSample(queryResults[i].VID),
                                m_index->GetSample(selections[j].headID));

                                // LOG(Helper::LogLevel::LL_Info,  "NNDist: %f Original: %f\n", nnDist, queryResults[i].Score);
                                if (nnDist <= queryResults[i].Dist)
                                {
                                    rngAccpeted = false;
                                    break;
                                }
                            }

                            if (!rngAccpeted)
                            {
                                continue;
                            }

                            selections[replicaCount].headID = queryResults[i].VID;
                            selections[replicaCount].fullID = VID;
                            selections[replicaCount].distance = queryResults[i].Dist;
							selections[replicaCount].order = (char)replicaCount;
                            ++replicaCount;
						}
						/*
						for (int i = 0; i < replicaCount; ++i)
						{
							LOG(Helper::LogLevel::LL_Info, "Insert headID: %d : Insert VID: %d\n", selections[i].headID, selections[i].fullID);
						}
						*/
						std::string postingList;
						//LOG(Helper::LogLevel::LL_Info, "Insert ID: %d, replica Num: %d\n", VID, replicaCount);
						for (int i = 0; i < replicaCount; ++i)
						{
							postingList.resize(0);
							postingList.clear();
							//LOG(Helper::LogLevel::LL_Info, "Insert VID: %d : Head Vector: %d\n", VID, selections[i].headID);
							db->Get(ReadOptions(), Helper::Serialize<int>(&selections[i].headID, 1), &postingList);
							int vectorNum = postingList.size()/ (sizeof(int) + m_vectorSize);
							std::shared_ptr<uint8_t> postingBuffer;
							postingBuffer.reset(new uint8_t[postingList.size() + m_vectorSize + sizeof(int)], std::default_delete<uint8_t[]>());
							uint8_t* bufferVoidPtr = reinterpret_cast<uint8_t*>(postingBuffer.get());
							memcpy(bufferVoidPtr, postingList.data(), postingList.size());
							bool covered = false;
							/*
							LOG(Helper::LogLevel::LL_Info, "Before Scan:\n");
							for (int j = 0; j < vectorNum; j++)
							{
								uint8_t* vectorInfo = postingBuffer.get() + j * (m_vectorSize + sizeof(int));
								int vectorID = *(reinterpret_cast<int*>(vectorInfo));
								LOG(Helper::LogLevel::LL_Info, " %d", vectorID);
							}
							LOG(Helper::LogLevel::LL_Info, "\n");
							*/
							for (int j = 0; j < vectorNum; j++)
							{
								uint8_t* vectorInfo = postingBuffer.get() + j * (m_vectorSize + sizeof(int));
								int vectorID = *(reinterpret_cast<int*>(vectorInfo));
								if (m_deletedID.Contains(vectorID)) 
								{
									//LOG(Helper::LogLevel::LL_Info, "Found Deleted ID: %d, replaced with ID: %d\n", vectorID, VID);
									covered = true;
									memcpy(vectorInfo, &VID, sizeof(int));
									memcpy(vectorInfo+sizeof(int), (uint8_t*)p_queryResults.GetTarget(), m_vectorSize);
									break;
								}
							}
							if (covered)
							{
								/*
								LOG(Helper::LogLevel::LL_Info, "After Insert:\n");
								for (int j = 0; j < vectorNum; j++)
								{
									uint8_t* vectorInfo = postingBuffer.get() + j * (m_vectorSize + sizeof(int));
									int vectorID = *(reinterpret_cast<int*>(vectorInfo));
									LOG(Helper::LogLevel::LL_Info, " %d", vectorID);
								}
								LOG(Helper::LogLevel::LL_Info, "\n");
								*/
								int datasize = postingList.size();
								bufferVoidPtr = reinterpret_cast<uint8_t*>(postingBuffer.get());
								postingList.resize(0);
								postingList.clear();
								postingList += Helper::Serialize<uint8_t>(bufferVoidPtr, datasize);
								db->Put(WriteOptions(), Helper::Serialize<int>(&selections[i].headID, 1), postingList);
								continue;
							}
							//LOG(Helper::LogLevel::LL_Info, "No Found Deleted ID, Insert ID: %d\n", VID);
							postingList += Helper::Serialize<int>(&VID, 1);
							postingList += Helper::Serialize<ValueType>(p_queryResults.GetTarget(), COMMON_OPTS.m_dim);
							vectorNum += 1;
							if (postingList.size() > m_postingPageLimit * p_pageSize)
							{
								//need to split and insert into headinedex
								//LOG(Helper::LogLevel::LL_Info, "PostingList Oversize, Need to Split\n");
								//extract out vector and id from postingList
								COMMON::Dataset<ValueType> smallSample;
								std::shared_ptr<uint8_t> vectorBuffer;
								//smallSample[i] -> VID
								std::vector<int> localindicesInsert;
								//localindices for kmeans smallSample[i] -> localindices[j] = i
								std::vector<int> localindices;
								localindicesInsert.resize(vectorNum);
								localindices.resize(vectorNum);
								bufferVoidPtr = reinterpret_cast<uint8_t*>(postingBuffer.get());
								memcpy(bufferVoidPtr, postingList.data(), postingList.size());
								vectorBuffer.reset(new uint8_t[m_vectorSize * vectorNum], std::default_delete<uint8_t[]>());
								bufferVoidPtr = reinterpret_cast<uint8_t*>(vectorBuffer.get());
								for (int j = 0; j < vectorNum; j++)
								{
									uint8_t* vectorInfo = postingBuffer.get() + j * (m_vectorSize + sizeof(int));
									localindicesInsert[j] = *(reinterpret_cast<int*>(vectorInfo));
									localindices[j] = j;
									memcpy(bufferVoidPtr, vectorInfo + sizeof(int), m_vectorSize);
									bufferVoidPtr += m_vectorSize;
								}
								smallSample.Initialize(vectorNum, COMMON_OPTS.m_dim, reinterpret_cast<ValueType*>(vectorBuffer.get()), false);
								//LOG(Helper::LogLevel::LL_Info, "Headid: %d Sample Vector Num: %d, Real Vector Num: %d\n", selections[i].headID, smallSample.R(), vectorNum);
								//k = 2, maybe we can change the split number
								SPTAG::COMMON::KmeansArgs<ValueType> args(m_k, smallSample.C(), (SizeType)localindicesInsert.size(), 1, m_index->GetDistCalcMethod());
								
								//int fBalanceFactor = SPTAG::COMMON::DynamicFactorSelect(smallSample, localindices, 0, (SizeType)localindices.size(), args);

								std::random_shuffle(localindices.begin(), localindices.end());

								int numClusters = SPTAG::COMMON::KmeansClustering(smallSample, localindices, 0, (SizeType)localindices.size(), args);

								if(numClusters <= 1)
								{
									//LOG(Helper::LogLevel::LL_Error, "Insert Stage:Very Close, First We Ignore spliting\n");
									db->Put(WriteOptions(), Helper::Serialize<int>(&selections[i].headID, 1), postingList);
									return ErrorCode::Success;
								}

								postingList.resize(0);
								postingList.clear();
								long long newHeadVID = -1;
								int first = 0;
								bool isHeadUpdate = false;
								std::vector<SizeType> fatherNodes;
								fatherNodes.clear();
								fatherNodes.emplace_back(selections[i].headID);
								for (int k = 0; k < m_k; k++) 
								{
									if (args.counts[k] == 0) continue;
									newHeadVID = -1;
									for (int j = 0; j < args.counts[k]; j++)
									{
										if (!isHeadUpdate && localindicesInsert[localindices[first + j]] == *m_vectorTranslateMap[selections[i].headID])
										{
											newHeadVID = selections[i].headID;
											isHeadUpdate = true;
										}
										postingList += Helper::Serialize<int>(&localindicesInsert[localindices[first + j]], 1);
                        				postingList += Helper::Serialize<ValueType>(smallSample[localindices[first + j]], COMMON_OPTS.m_dim);
									}
									if (newHeadVID == -1)
									{
										if (!isHeadUpdate && k == m_k-1)
										{
											//the last new head vector and the privous father head vector is deleted
											//add this postinglist into head vector
											newHeadVID = selections[i].headID;
										}
										else
										{
											//LOG(Helper::LogLevel::LL_Info, "Insert new head vector\n");
											newHeadVID = localindicesInsert[localindices[first + args.counts[k] - 1]];
											//BUG: newHeadVID maybe a exist head vector
											m_vectorTranslateMap.AddBatch(&newHeadVID, 1);
											newHeadVID = m_vectorTranslateMap.R() - 1;
											m_split_num++;
											m_index->AddHeadIndex(smallSample[localindices[first + args.counts[k] - 1]], 1, COMMON_OPTS.m_dim, fatherNodes);
										}
									}
									LOG(Helper::LogLevel::LL_Info, "Headid: %d split into : %d\n", selections[i].headID, newHeadVID);
									db->Put(WriteOptions(), Helper::Serialize<int>(&newHeadVID, 1), postingList);
									postingList.resize(0);
									postingList.clear();
									first += args.counts[k];
								}
							} else
							{
								db->Put(WriteOptions(), Helper::Serialize<int>(&selections[i].headID, 1), postingList);
							}
						}
					} else {
						LOG(Helper::LogLevel::LL_Error, "Only Support BKT Update");
						return ErrorCode::Undefined;
					}
					return ErrorCode::Success;
				}

				ErrorCode Insert(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats)
				{
					m_deletedID.AddBatch(1);
					m_index->SearchIndex(p_queryResults);
					auto ret = InsertPostingList(p_queryResults, p_stats, m_vectornum);
					m_vectornum++;
					return ret;
				}

				ErrorCode Delete(const SizeType& p_id) {
            		std::shared_lock<std::shared_timed_mutex> sharedlock(m_dataDeleteLock);
            		if (m_deletedID.Insert(p_id)) return ErrorCode::Success;
            		return ErrorCode::VectorNotFound;
        		}

        		ErrorCode Delete(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats) 
				{
					p_queryResults.GetTarget();
					Search(p_queryResults, p_stats);
#pragma omp parallel for schedule(dynamic)
            		for (SizeType i = 0; i < p_queryResults.GetResultNum(); i++) {
                    	if (p_queryResults.GetResult(i)->Dist < 1e-6) {
                        	Delete(p_queryResults.GetResult(i)->VID);
                    	}
            		}
            		return ErrorCode::Success;
        		}

				ErrorCode Save()
				{
					m_deletedID.Save(COMMON_OPTS.m_deleteID);
					m_index->SaveIndex(COMMON_OPTS.m_headIndexFolder);
					return ErrorCode::Success;
				}

				ErrorCode CalDBDist(std::string& storefile)
				{
					std::string postingList;
					auto ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize(storefile.c_str(), std::ios::binary | std::ios::out))
                    {
                        LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", storefile.c_str());
                        exit(1);
                    }
					for (int i = 0; i < m_vectorTranslateMap.R(); i++)
					{
						postingList.clear();
						db->Get(ReadOptions(), Helper::Serialize<int>(&i, 1), &postingList);
						int vectorNum = postingList.size()/ (sizeof(int) + m_vectorSize);
						std::shared_ptr<uint8_t> postingBuffer;
						postingBuffer.reset(new uint8_t[postingList.size() + m_vectorSize + sizeof(int)], std::default_delete<uint8_t[]>());
						uint8_t* bufferVoidPtr = reinterpret_cast<uint8_t*>(postingBuffer.get());
						memcpy(bufferVoidPtr, postingList.data(), postingList.size());
						int size = 0;
						for (int j = 0; j < vectorNum; j++)
						{
							uint8_t* vectorInfo = postingBuffer.get() + j * (m_vectorSize + sizeof(int));
							int vectorID = *(reinterpret_cast<int*>(vectorInfo));
							if (m_deletedID.Contains(vectorID)) 
							{
								continue;
							}
							size++;
						}
						if (ptr->WriteBinary(sizeof(int), (char*)&i) != sizeof(int)) {
                        	LOG(Helper::LogLevel::LL_Error, "Fail to write head");
                        	exit(1);
                    	}
						if (ptr->WriteBinary(sizeof(int), (char*)&vectorNum) != sizeof(int)) {
                        	LOG(Helper::LogLevel::LL_Error, "Fail to write posting size");
                        	exit(1);
                    	}
						if (ptr->WriteBinary(sizeof(int), (char*)&size) != sizeof(int)) {
                        	LOG(Helper::LogLevel::LL_Error, "Fail to write real posting size");
                        	exit(1);
                    	}
					}
					return ErrorCode::Success;
				}

				void Search(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats)
				{
					//LARGE_INTEGER qpcStartTime;
					//LARGE_INTEGER qpcEndTime;
					//QueryPerformanceCounter(&qpcStartTime);
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
						//LOG(Helper::LogLevel::LL_Info, "Adding PostingList\n");
						for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
						{
							auto res = p_queryResults.GetResult(i);
							if (res->VID != -1)
							{
								auto_ws->m_postingIDs.emplace_back(res->VID);
								p_stats.m_headAndDist[res->VID] = res->Dist;
							}
						}
						const uint32_t postingListCount = static_cast<uint32_t>(auto_ws->m_postingIDs.size());
                   		m_extraSearcher->InitWorkSpace(auto_ws, postingListCount);
					}

					if (!COMMON_OPTS.m_addHeadToPost && m_vectorTranslateMap.Name() == "Head")
					{
						//LOG(Helper::LogLevel::LL_Info, "Translate headvector\n");
						for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
						{
							auto res = p_queryResults.GetResult(i);
							if (res->VID != -1)
							{
								//LOG(Helper::LogLevel::LL_Info, "head vector previous id: %d\n", res->VID);
								res->VID = static_cast<int>(*m_vectorTranslateMap[res->VID]);
								//LOG(Helper::LogLevel::LL_Info, "head vector normal id: %d\n", res->VID);
								if (auto_ws->m_deduper.CheckAndSet(res->VID) || m_deletedID.Contains(res->VID))
								{
									//LOG(Helper::LogLevel::LL_Info, "contain\n");
									res->VID = -1;
									res->Dist = SPTAG::MaxDist;
									res->Meta.Clear();
								}
								//LOG(Helper::LogLevel::LL_Info, "get next head vector\n");
							}
						}
					} else if (COMMON_OPTS.m_addHeadToPost && m_clearHead)
					{
						//the head vector will be count at extraSearcher
						p_queryResults.Reset();
					}

					if (nullptr != m_extraSearcher)
					{
						p_queryResults.Reverse();
						//LOG(Helper::LogLevel::LL_Info, "Into PostingList Search\n");
						m_extraSearcher->Search(auto_ws, p_queryResults, m_index, p_stats, m_deletedID);
						RetWs(auto_ws);
					}

					ExEndingTime = sw.getElapsedMs();
					//QueryPerformanceCounter(&qpcEndTime);
					//unsigned __int32 latency = static_cast<unsigned __int32>(static_cast<double>(qpcEndTime.QuadPart - qpcStartTime.QuadPart) * 1000000 / g_systemPerfFreq.QuadPart);

					p_stats.m_exLatency = ExEndingTime - EndingTime;
					p_stats.m_totalSearchLatency = ExEndingTime - StartingTime;
					p_stats.m_totalLatency = p_stats.m_totalSearchLatency;
					//p_stats.m_totalLatency = latency * 1.0;
				}

				void Search4ANNIndexTestTool(COMMON::QueryResultSet<ValueType>& p_queryResults)
				{
					m_index->SearchIndex(p_queryResults);

					ExtraWorkSpace* auto_ws = nullptr;
					auto_ws = GetWs();
					auto_ws->m_postingIDs.clear();

					for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
					{
						auto res = p_queryResults.GetResult(i);
						if (res->VID != -1)
						{
							auto_ws->m_postingIDs.emplace_back(res->VID);
							res->VID = static_cast<int>(*m_vectorTranslateMap[res->VID]);
						}
					}

					p_queryResults.Reverse();

					m_extraSearcher->Search(auto_ws, p_queryResults, m_index);
					RetWs(auto_ws);
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

				void SearchAsync(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback)
				{
					p_stats.m_searchRequestTime = std::chrono::steady_clock::now();

					SearchAsyncJob* curJob = new SearchAsyncJob(this, p_queryResults, p_stats, p_callback);

					m_threadPool->add(curJob);
				}


				void SetHint(int p_threadNum, int p_resultNum, bool p_asyncCall, const Options& p_opts)
				{
					LOG(Helper::LogLevel::LL_Info, "ThreadNum: %d, ResultNum: %d, AsyncCall: %d\n", p_threadNum, p_resultNum, p_asyncCall ? 1 : 0);

					m_internalResultNum = p_resultNum;

					m_replicaCount = p_opts.m_replicaCount;

					m_postingPageLimit = p_opts.m_postingPageLimit;

					m_vectorSize = COMMON_OPTS.m_dim * sizeof(ValueType);

					m_split_num = 0;

					if (p_asyncCall)
					{
						m_threadPool.reset(new Helper::ThreadPool());
						m_threadPool->init(p_threadNum);
					}
				}

				void setSplitZero()
				{
					m_split_num = 0;
				}

				int getSplitNum()
				{
					return m_split_num;
				}

				std::shared_ptr<VectorIndex> HeadIndex() {
					return m_index;
				}

				ExtraWorkSpace* GetWs() {
					ExtraWorkSpace* ws = nullptr;
					if (!m_workspaces.pop(ws)) {
						ws = new ExtraWorkSpace();
					}

					return ws;
				}

				void RetWs(ExtraWorkSpace* ws) {
					if (ws != nullptr)
					{
						m_workspaces.push(ws);
					}
				}

			protected:
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

				std::shared_ptr<VectorIndex> m_index;

				SPTAG::COMMON::Dataset<long long> m_vectorTranslateMap;
				//std::unique_ptr<long long[]> m_vectorTranslateMap;
				//std::vector<long long> m_vectorTranslateMap;

				std::unique_ptr<IExtraSearcher<ValueType>> m_extraSearcher;

				std::unique_ptr<Helper::ThreadPool> m_threadPool;

				std::atomic<std::int32_t> m_tids;

				//for real-time service and release
				COMMON::Labelset m_deletedID;

				//delete mutex lock
				std::shared_timed_mutex m_dataDeleteLock;

				int m_internalResultNum;

				int m_replicaCount;

				int m_postingPageLimit;

				int m_vectorSize;

				int m_k = 2;

				int m_vectornum;

				bool m_clearHead = true;

				int m_split_num;

				boost::lockfree::stack<ExtraWorkSpace*> m_workspaces;
			};
		}
	}
}
