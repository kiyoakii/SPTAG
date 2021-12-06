#pragma once

#ifdef _MSC_VER
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcher.h"
#else // non windows
#include "inc/SSDServing/VectorSearch/IExtraSearcherLinux.h"
#include "inc/SSDServing/VectorSearch/ExtraFullGraphSearcherLinux.h"
#endif

#include "inc/Core/BKT/Index.h"
#include "inc/Core/Common/BKTree.h"
#include "inc/Core/Common/Dataset.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Core/VectorIndex.h"
#include "inc/Helper/ThreadPool.h"
#include "inc/Helper/Concurrent.h"
#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"
#include "inc/SSDServing/IndexBuildManager/Utils.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"
#include "inc/SSDServing/VectorSearch/PersistentBuffer.h"

#include <boost/lockfree/stack.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/lockfree/policies.hpp>

#include <atomic>
#include <shared_mutex>
#include <chrono>
#include <thread>

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

			enum taskType { del=0, ins=1 };

			struct Task {
				Task() : type(ins), id(0), appendNum(0), part(nullptr) {}
				Task(taskType type, SizeType id) : type(type), id(id), appendNum(0), part(nullptr) {}
				Task(taskType type, SizeType id, uint32_t num, std::string* part) : 
					type(type), id(id), appendNum(num), part(part) {}
				Task(const Task &t)
				{
					type = t.type;
					id = t.id;
					appendNum = t.appendNum;
					part = t.part;
				}
				taskType    type;
				SizeType   	id;
				uint32_t    appendNum;
				std::string* part;
			};

			template <typename ValueType>
			class SearchDefault
			{
			public:
				SearchDefault()
					: m_workspaces(128), currentTasks(60000), m_skipped(0),
					running(new uint8_t[1000000000]), splitRoute(new std::pair<SizeType, SizeType>[1000000000])
				{
					m_tids = 0;
					m_replicaCount = 4;
					//QueryPerformanceFrequency(&g_systemPerfFreq);
					for (int i = 0; i < 1000000000; i++) {
						running[i] = 0;
						splitRoute[i] = std::make_pair<SizeType, SizeType>(0, 0);
					}
				}

				~SearchDefault()
				{
					ExtraWorkSpace* context;
					while (m_workspaces.pop(context))
					{
						delete context;
					}
					delete running;
					delete splitRoute;
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

					int vectornum;

                    input.read(reinterpret_cast<char*>(&vectornum), sizeof(vectornum));

					m_vectornum.store(vectornum);

					LOG(Helper::LogLevel::LL_Info, "Current vector num: %d.\n", m_vectornum.load());

					input.read(reinterpret_cast<char*>(&m_posting_num), sizeof(m_posting_num));

					LOG(Helper::LogLevel::LL_Info, "Current posting num: %d.\n", m_posting_num);

					m_postingSizes.resize(m_posting_num);
					m_postingRadius.resize(m_posting_num);

					input.read(reinterpret_cast<char*>(m_postingSizes.data()), sizeof(int) * m_posting_num);
					input.read(reinterpret_cast<char*>(m_postingRadius.data()), sizeof(float) * m_posting_num);

					input.close();

					m_vectornum_last = m_vectornum.load();
					
					calAvgPostingSize();

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
						if (p_config.m_persistentBufferPath.empty())
						{
							LOG(Helper::LogLevel::LL_Error, "Error in read persistent buffer path.\n");
						}
						m_persistentBuffer.reset(new PersistentBuffer(p_config.m_persistentBufferPath));
					}
					m_clearHead = !p_config.m_buildSsdIndex;
				}

				void testPrint(int headID)
				{
					std::string postingList;
					//LOG(Helper::LogLevel::LL_Info, "reading posting\n");
					db->Get(ReadOptions(), Helper::Serialize<int>(&headID, 1), &postingList);
					//LOG(Helper::LogLevel::LL_Info, "reading posting finish\n");
					int vectorNum = postingList.size() / (sizeof(int) + m_vectorSize);
					uint8_t* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
					for (int i = 0; i < vectorNum; i++)
					{
						uint8_t* vectorId = postingP + i * (sizeof(int) + m_vectorSize);
						LOG(Helper::LogLevel::LL_Info, "vector ID: %d ", *(reinterpret_cast<int*>(vectorId)));
					}
					LOG(Helper::LogLevel::LL_Info, "\n");
				}

				void Rebuild()
				{
					// m_index should have a virtual rebuild method, now I just cast
					LOG(Helper::LogLevel::LL_Info, "Rebuild Head Index...\n");
					auto bkt = dynamic_cast<SPTAG::BKT::Index<ValueType>*>(m_index.get());
					bkt->Rebuild();
					LOG(Helper::LogLevel::LL_Info, "Finish Rebuild Head Index...\n");
				}

				void setPair(std::pair<SizeType, SizeType>& p, int k, SizeType id) {
					if (k == 1) {
						p.first = id;
					} else if (k == 2) {
						p.second = id;
					}
				}

				ErrorCode Split(const SizeType headID, int appendNum, std::string& appendPosting)
				{
					//int father = m_index->GetFatherID(headID);
					//LOG(Helper::LogLevel::LL_Info, "headID: %d, fatherID: %d\n", headID, father);
					std::string postingList;
					db->Get(ReadOptions(), Helper::Serialize<int>(&headID, 1), &postingList);
					postingList += appendPosting;
					uint8_t* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
					int postVectorNum =  postingList.size() / (m_vectorSize + sizeof(int));
					// reinterpret postingList to vectors and IDs
					COMMON::Dataset<ValueType> smallSample;  // smallSample[i] -> VID
					std::shared_ptr<uint8_t> vectorBuffer(new uint8_t[m_vectorSize * postVectorNum], std::default_delete<uint8_t[]>());
					std::vector<int> localindicesInsert(postVectorNum);  // smallSample[i] = j <-> localindices[j] = i
					std::vector<int> localindices(postVectorNum);
					auto vectorBuf = vectorBuffer.get();
					int realVectorNum = postVectorNum;
					int index = 0;
					//LOG(Helper::LogLevel::LL_Info, "Scanning\n");
					for (int j = 0; j < postVectorNum; j++)
					{
						uint8_t* vectorId = postingP + j * (m_vectorSize + sizeof(int));
						//LOG(Helper::LogLevel::LL_Info, "vector index/total:id: %d/%d:%d\n", j, m_postingSizes[selections[i].headID], *(reinterpret_cast<int*>(vectorId)));
						if (m_deletedID.Contains(*(reinterpret_cast<int*>(vectorId)))) {
							realVectorNum--;
						} else {
							localindicesInsert[index] = *(reinterpret_cast<int*>(vectorId));
							localindices[index] = index;
							index++;
							memcpy(vectorBuf, vectorId + sizeof(int), m_vectorSize);
							vectorBuf += m_vectorSize;
						}
					}
					if (realVectorNum < postVectorNum * 0.9)
					{
						postingList.clear();
						float r = 0.f;
						for (int j = 0; j < realVectorNum; j++)
						{
							postingList += Helper::Serialize<int>(&localindicesInsert[j], 1);
							postingList += Helper::Serialize<ValueType>(vectorBuffer.get() + j * m_vectorSize, COMMON_OPTS.m_dim);
							auto dist = m_index->ComputeDistance(vectorBuffer.get() + j * m_vectorSize, m_index->GetSample(headID));
							r = std::max<float>(r, dist);
						}
						m_postingSizes[headID] = realVectorNum;
						m_postingRadius[headID] = r;
						db->Put(WriteOptions(), Helper::Serialize<int>(&headID, 1), postingList);
						return ErrorCode::Success;
					}
					//LOG(Helper::LogLevel::LL_Info, "Resize\n");
					localindicesInsert.resize(realVectorNum);
					localindices.resize(realVectorNum);
					smallSample.Initialize(realVectorNum, COMMON_OPTS.m_dim, reinterpret_cast<ValueType*>(vectorBuffer.get()), false);
					
					//LOG(Helper::LogLevel::LL_Info, "Headid: %d Sample Vector Num: %d, Real Vector Num: %d\n", selections[i].headID, smallSample.R(), realVectorNum);
					
					// k = 2, maybe we can change the split number
					SPTAG::COMMON::KmeansArgs<ValueType> args(m_k, smallSample.C(), (SizeType)localindicesInsert.size(), 1, m_index->GetDistCalcMethod());
					std::random_shuffle(localindices.begin(), localindices.end());
					int numClusters = SPTAG::COMMON::KmeansClustering(smallSample, localindices, 0, (SizeType)localindices.size(), args);
					if (numClusters <= 1)
					{
						postingList.clear();
						// float r = 0.f;
						for (int j = 0; j < realVectorNum; j++)
						{
							postingList += Helper::Serialize<int>(&localindicesInsert[j], 1);
							postingList += Helper::Serialize<ValueType>(vectorBuffer.get() + j * m_vectorSize, COMMON_OPTS.m_dim);
							// auto dist = m_index->ComputeDistance(vectorBuffer.get() + j * m_vectorSize, m_index->GetSample(headID));
							// r = std::max<float>(r, dist);
						}
						m_postingSizes[headID] = realVectorNum;
						m_postingRadius[headID] = args.clusterDist[0];
						db->Put(WriteOptions(), Helper::Serialize<int>(&headID, 1), postingList);
						return ErrorCode::Success;
					}
					
					long long newHeadVID = -1;
					int first = 0;
					//std::vector<SizeType> fatherNodes;
					//fatherNodes.emplace_back(father);
					std::pair<SizeType, SizeType> p;
					for (int k = 0; k < m_k; k++) 
					{
						int begin, end = 0;
						std::string postingList;
						if (args.counts[k] == 0)	continue;

						// LOG(Helper::LogLevel::LL_Info, "Insert new head vector\n");
						// Notice: newHeadVID maybe a exist head vector

						m_index->AddHeadIndexId(smallSample[args.clusterIdx[k]], 1, COMMON_OPTS.m_dim, &begin, &end);
						newHeadVID = begin;
						setPair(p, k+1, newHeadVID);

						//LOG(Helper::LogLevel::LL_Info, "Headid: %d split into : %d\n", headID, newHeadVID);
						float MaxClusterDist = args.clusterDist[k];
						for (int j = 0; j < args.counts[k]; j++)
						{
							postingList += Helper::Serialize<int>(&localindicesInsert[localindices[first + j]], 1);
							postingList += Helper::Serialize<ValueType>(smallSample[localindices[first + j]], COMMON_OPTS.m_dim);
							float newDist = m_index->ComputeDistance(smallSample[args.clusterIdx[k]], smallSample[localindices[first + j]]);
							if (newDist > MaxClusterDist)
							{
								MaxClusterDist = newDist;
							}
						}
						db->Put(WriteOptions(), Helper::Serialize<int>(&newHeadVID, 1), postingList);
						first += args.counts[k];

						{
							std::lock_guard<std::mutex> lock(m_headAddLock);
							/*
							std::unique_ptr<std::atomic_int> newAtomicSize(new std::atomic_int(args.counts[k]));
							m_postingSizes.push_back(std::move(newAtomicSize));
							*/
							m_postingSizes.resize(m_index->GetNumSamples());
							m_postingRadius.resize(m_index->GetNumSamples());
							m_postingSizes[newHeadVID] = args.counts[k];
							m_postingRadius[newHeadVID] = MaxClusterDist;
							m_posting_num++;
							m_split_num++;
						}
						m_index->AddHeadIndexIdx(begin, end);
					}
					splitRoute[headID] = p;

					// delete from BKT and RNG
					auto bktIndex = dynamic_cast<SPTAG::BKT::Index<ValueType>*>(m_index.get());
					bktIndex->DeleteIndex(headID);
					// delete from disk
					db->Delete(WriteOptions(), Helper::Serialize<int>(&headID, 1));
					//delete m_postingSizes[headID].get();
					//m_postingSizes[headID].reset(nullptr);
					return ErrorCode::Success;
				}

				ErrorCode Append(const SizeType headID, int appendNum, std::string* appendPosting)
				{
					if (m_postingSizes[headID] + appendNum > m_postingVectorLimit){
						Split(headID, appendNum, *appendPosting);
					} else {
						//LOG(Helper::LogLevel::LL_Info, "Merge: headID: %d, appendNum:%d\n", headID, appendNum);
						db->Merge(WriteOptions(), Helper::Serialize<int>(&headID, 1), *appendPosting);
						//LOG(Helper::LogLevel::LL_Info, "here\n");
						m_postingSizes[headID] += appendNum;
						//LOG(Helper::LogLevel::LL_Info, "Scan\n");
						uint8_t* postingP = reinterpret_cast<uint8_t*>(&appendPosting->front()) + sizeof(int);
						float r = m_postingRadius[headID];
						for (int i = 0; i < appendNum; i++) {
							r = std::max<float>(r, m_index->ComputeDistance(m_index->GetSample(headID), postingP));
							postingP += m_vectorSize + sizeof(int);
						}
						m_postingRadius[headID] = r;
					}
					running[headID] = 0;
					//LOG(Helper::LogLevel::LL_Info, "Delete\n");
					delete appendPosting;
					//LOG(Helper::LogLevel::LL_Info, "Finish\n");

					return ErrorCode::Success;
				}

				ErrorCode Delete(const SizeType& p_id) {
            		std::shared_lock<std::shared_timed_mutex> sharedlock(m_dataDeleteLock);
            		if (m_deletedID.Insert(p_id)) return ErrorCode::Success;
            		return ErrorCode::VectorNotFound;
        		}

        		ErrorCode Delete(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats) 
				{
					Search(p_queryResults, p_stats);
#pragma omp parallel for schedule(dynamic)
            		for (SizeType i = 0; i < p_queryResults.GetResultNum(); i++) {
                    	if (p_queryResults.GetResult(i)->Dist < 1e-6) {
                        	Delete(p_queryResults.GetResult(i)->VID);
                    	}
            		}
            		return ErrorCode::Success;
        		}

				SizeType Updater(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, int* VID)
				{
					*VID = m_vectornum.fetch_add(1);
					{
						std::lock_guard<std::mutex> lock(m_dataAddLock);
						m_deletedID.AddBatch(1);
					}
					m_index->SearchIndex(p_queryResults);
					if (COMMON_OPTS.m_indexAlgoType != IndexAlgoType::BKT) {
						LOG(Helper::LogLevel::LL_Error, "Only Support BKT Update");
						return -1;
					}
					int replicaCount = 0;
					BasicResult* queryResults = p_queryResults.GetResults();
					std::vector<EdgeInsert> selections(static_cast<size_t>(m_replicaCount));
					for (int i = 0; i < p_queryResults.GetResultNum() && replicaCount < m_replicaCount; ++i)
					{
						if (queryResults[i].VID == -1) {
							break;
						}
						// LOG(Helper::LogLevel::LL_Info, "Head Vector: %d\n", queryResults[i].VID);

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
							continue;
						selections[replicaCount].headID = queryResults[i].VID;
						selections[replicaCount].fullID = *VID;
						selections[replicaCount].distance = queryResults[i].Dist;
						selections[replicaCount].order = (char)replicaCount;
						++replicaCount;
					}
					char insertCode = 0;
					SizeType assignID = 0;
					for (int i = 0; i < replicaCount; i++)
					{
						std::string assignment;
						assignment += Helper::Serialize<char>(&insertCode, 1);
						assignment += Helper::Serialize<int>(&selections[i].headID, 1);
						assignment += Helper::Serialize<int>(VID, 1);
						assignment += Helper::Serialize<ValueType>(p_queryResults.GetTarget(), COMMON_OPTS.m_dim);
						assignID = m_persistentBuffer->PutAssignment(assignment);
					}
					return assignID;
				}

				//delete updater
				ErrorCode Updater(const SizeType& p_id)
				{
					char deleteCode = 1;
					int VID = p_id;
					std::string assignment;
					assignment += Helper::Serialize<char>(&deleteCode, 1);
					assignment += Helper::Serialize<int>(&VID, 1);
					m_persistentBuffer->PutAssignment(assignment);
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
						if (m_index->ContainSample(i)) {
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
					}
					return ErrorCode::Success;
				}

				void Search_VecLimit(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats)
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
						int totalVectors = 0;
						for (int i = 0; i < p_queryResults.GetResultNum(); ++i)
						{
							auto res = p_queryResults.GetResult(i);
							if (res->VID != -1)
							{
								auto_ws->m_postingIDs.emplace_back(res->VID);
								totalVectors += m_postingSizes[res->VID];
								if (totalVectors > m_searchVectorLimit)
								{
									break;
								}
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
						int totalVectors = 0;
						float currentR = p_queryResults.GetResult(m_resultNum - 1)->Dist;
						// First, add topK (resultNum) headVectors.
						// TopK actually represents vector number to be returned,
						// but here it acts as a bound to form a search area
						// LOG(Helper::LogLevel::LL_Info, "m_resultNum = %d\n", m_resultNum);
						for (int i = 0; i < m_resultNum && totalVectors <= m_searchVectorLimit; ++i) {
							auto res = p_queryResults.GetResult(i);
							if (res->VID != -1) {
								auto_ws->m_postingIDs.emplace_back(res->VID);
								currentR = res->Dist;
								totalVectors += m_postingSizes[res->VID];
							}
						}
						for (int i = m_resultNum; i < p_queryResults.GetResultNum() && totalVectors <= m_searchVectorLimit; ++i) {
							auto res = p_queryResults.GetResult(i);
							if (res->VID != -1 && res->Dist - m_postingRadius[res->VID] <= currentR) {
								auto_ws->m_postingIDs.emplace_back(res->VID);
								totalVectors += m_postingSizes[res->VID];
							} else {
								m_skipped++;
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

				class AppendAsyncJob : public SPTAG::Helper::ThreadPool::Job
				{
				private:
					SearchDefault* m_processor;
					const SizeType headID;
					int appendNum;
					std::string* appendPosting;
					std::function<void()> m_callback;
				public:
					AppendAsyncJob(SearchDefault* p_processor,
						const SizeType headID, int appendNum, std::string* appendPosting, std::function<void()> p_callback)
						: m_processor(p_processor),
						headID(headID), appendNum(appendNum), appendPosting(appendPosting), m_callback(p_callback) {}

					~AppendAsyncJob() {}

					void exec() {
						m_processor->ProcessAsyncAppend(headID, appendNum, appendPosting, std::move(m_callback));
					}
				};

				void AppendAsync(const SizeType headID, int appendNum, std::string* appendPosting, std::function<void()> p_callback=nullptr)
				{
					AppendAsyncJob* curJob = new AppendAsyncJob(this, headID, appendNum, appendPosting, p_callback);
					m_appendThreadPool->add(curJob);
				}

				void PushAsync(const Task t)
				{
					PushAsyncJob* curJob = new PushAsyncJob(this, t);
					m_appendThreadPool->add(curJob);
				}

				class DeleteAsyncJob : public SPTAG::Helper::ThreadPool::Job
				{
				private:
					SearchDefault* m_processor;
					const SizeType& p_id;
					std::function<void()> m_callback;
				public:
					DeleteAsyncJob(SearchDefault* p_processor,
						const SizeType& p_id, std::function<void()> p_callback)
						: m_processor(p_processor),
						p_id(p_id), m_callback(p_callback) {}

					~DeleteAsyncJob() {}

					void exec() {
						m_processor->ProcessAsyncDelete(p_id, std::move(m_callback));
					}
				};

				class PushAsyncJob : public SPTAG::Helper::ThreadPool::Job
				{
				private:
					SearchDefault* m_processor;
					Task t;
				public:
					PushAsyncJob(SearchDefault* p_processor, Task t)
						: m_processor(p_processor), t(t){}

					~PushAsyncJob() {}

					void exec() {
						m_processor->ProcessAsyncPush(t);
					}
				};

				void DeleteAsync(const SizeType& p_id, std::function<void()> p_callback=nullptr)
				{
					DeleteAsyncJob* curJob = new DeleteAsyncJob(this, p_id, p_callback);

					m_deleteThreadPool->add(curJob);
				}


				void SetHint(int p_threadNum, int p_resultNum, bool p_asyncCall, const Options& p_opts)
				{
					LOG(Helper::LogLevel::LL_Info, "ThreadNum: %d, ResultNum: %d, AsyncCall: %d\n", p_threadNum, p_resultNum, p_asyncCall ? 1 : 0);

					m_internalResultNum = p_resultNum;

					m_resultNum = p_opts.m_resultNum;

					m_replicaCount = p_opts.m_replicaCount;

					m_postingPageLimit = p_opts.m_postingPageLimit;

					m_vectorSize = COMMON_OPTS.m_dim * sizeof(ValueType);

					m_postingVectorLimit = m_postingPageLimit * p_pageSize / (sizeof(int) + m_vectorSize);

					m_split_num = 0;

					m_k = p_opts.m_k;

					m_searchVectorLimit = p_opts.m_searchVectorLimit;

					if (p_asyncCall)
					{
						m_threadPool.reset(new Helper::ThreadPool());
						m_threadPool->init(p_threadNum);
					}
				}

				void updaterSetup()
				{
					appliedAssignment = 0;
					finishedAssignment = 0;
					m_appendThreadPool.reset(new Helper::ThreadPool());
					m_appendThreadPool->init(m_appendThreadNum);
					m_deleteThreadPool.reset(new Helper::ThreadPool());
					m_deleteThreadPool->init(m_deleteThreadNum);
					m_pushThreadPool.reset(new Helper::ThreadPool());
					m_pushThreadPool->init(m_pushThreadNum);
					m_dispatcher_running_flag.test_and_set();
					auto scanThread = std::thread(&SearchDefault::scanner, this);
					scanThread.detach();
					auto DispatchLoopThread = std::thread(&SearchDefault::dispatcher, this);
					DispatchLoopThread.detach();
				}

				void setSplitZero()
				{
					m_split_num = 0;
				}

				int getSplitNum()
				{
					return m_split_num;
				}

				int getVecNum()
				{
					return m_vectornum;
				}

				void calAvgPostingSize()
				{
					int total = 0;
					int deleted = 0;
					for (int i = 0; i < m_postingSizes.size(); i++)
					{
						if (m_index->ContainSample(i)) total += m_postingSizes[i];
						else deleted++;
					}
					m_postingSize_avg = total / (m_posting_num - deleted);
				}

				void setSearchLimit(int topK)
				{
					m_searchVectorLimit = topK * m_postingSize_avg;
					LOG(Helper::LogLevel::LL_Info, "Search Vector Limit: %d\n", m_searchVectorLimit);
				}

				std::shared_ptr<VectorIndex> HeadIndex() 
				{
					return m_index;
				}

				ExtraWorkSpace* GetWs() 
				{
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

				void scanner()
				{
					while (true) {
						bool noAssignment = true;
						int currentAssignmentID = m_persistentBuffer->GetCurrentAssignmentID();
						int scanNum = std::min<int>(appliedAssignment + 1000, currentAssignmentID);
						if (scanNum != appliedAssignment) {
							noAssignment = false;
						} else if (!m_dispatcher_running_flag.test_and_set()) {
							return;
						}

						std::map<SizeType, std::string*> newPart;
						std::set<SizeType> deletedVector;
						for (int i = appliedAssignment; i < scanNum; i++) {
							std::string assignment;
							m_persistentBuffer->GetAssignment(i, &assignment);
							uint8_t* postingP = reinterpret_cast<uint8_t*>(&assignment.front());
							char code = *(reinterpret_cast<char*>(postingP));
							if (code == 0) {
								// insert
								uint8_t* headPointer = postingP + sizeof(char);
								int32_t headID = *(reinterpret_cast<int*>(headPointer));
								if (newPart.find(headID) == newPart.end()) {
									newPart[headID] = new std::string(Helper::Serialize<uint8_t>(headPointer + sizeof(int), m_vectorSize + sizeof(int)));
								} else {
									*newPart[headID] += Helper::Serialize<uint8_t>(headPointer + sizeof(int), m_vectorSize + sizeof(int));
								}
							} else {
								// delete
								uint8_t* vectorPointer = postingP + sizeof(char);
								int VID = *(reinterpret_cast<int*>(vectorPointer));
								deletedVector.insert(VID);
							}
						}

						for (auto iter = newPart.begin(); iter != newPart.end(); iter++) {
							int appendNum = (*iter->second).size() / (m_vectorSize + sizeof(int));
							while(!currentTasks.push(Task(ins, iter->first, appendNum, iter->second))) {}
						}
						for (auto iter = deletedVector.begin(); iter != deletedVector.end(); iter++) {
							// while(!currentTasks.push(Task(del, *iter))) {}
							m_deletedID.Insert(*iter);
						}
						
						appliedAssignment = scanNum;
						if (noAssignment) {
							std::this_thread::sleep_for(std::chrono::milliseconds(100));
						} else {
							//LOG(Helper::LogLevel::LL_Info, "Process Append Assignments: %d, Delete Assignments: %d\n", newPart.size(), deletedVector.size());
						}
					}
				}

				std::vector<SizeType> traceSplitRoute(SizeType id) 
				{
					if (splitRoute[id] == std::make_pair<SizeType, SizeType>(0, 0)) {
						return std::vector<SizeType>{ id };
					} else {
						auto pair = splitRoute[id];
						std::vector<SizeType> ans;
						auto p1 = traceSplitRoute(pair.first);
						auto p2 = traceSplitRoute(pair.second);
						
						ans.insert(ans.end(), p1.begin(), p1.end());
						ans.insert(ans.end(), p2.begin(), p2.end());
						return ans;
					}
				}

				void dispatcher()
				{
					while (true) {
						Task task;
						while (!currentTasks.pop(task))
						{
							// no assignment sleep
							if (!m_dispatcher_running_flag.test_and_set()) {
								return;
							}
							std::this_thread::sleep_for(std::chrono::milliseconds(100));
						}
						// if (task.type == del) {
						// 	DeleteAsync(task.id);
						// } else {
						if (running[task.id] == 0) {
							running[task.id] = 1;
							if (!m_index->ContainSample(task.id)) {
								running[task.id] = 0;
								std::vector<SizeType> newIDs = traceSplitRoute(task.id);
								auto dis = [this, task](SizeType a) -> float { return m_index->ComputeDistance(m_index->GetSample(task.id), m_index->GetSample(a)); };
								auto sortFunc = [dis](SizeType a, SizeType b) -> bool { return dis(a) < dis(b); };
								std::sort(newIDs.begin(), newIDs.end(), sortFunc);
								if (running[newIDs[0]] == 1) {
									// new head is appending
									task.id = newIDs[0];
									if (task.id == 0) {
										LOG(Helper::LogLevel::LL_Info, "BUG: id = 0\n");
									}
									// if (!currentTasks.push(task)) {
									// 	LOG(Helper::LogLevel::LL_Error, "Lockfree queue capacity not enough!");
									// }
									PushAsync(task);
								} else {
									running[newIDs[0]] = 1;
									AppendAsync(newIDs[0], task.appendNum, task.part);
								}
							} else {
								AppendAsync(task.id, task.appendNum, task.part);
							}
						} else {
							// if (!currentTasks.push(task)) {
							// 	LOG(Helper::LogLevel::LL_Error, "Lockfree queue capacity not enough!");
							// }
							PushAsync(task);
						}
						// }
					}
				}

				bool checkAllTaskesIsFinish()
				{
				 	int currentAssignmentID = m_persistentBuffer->GetCurrentAssignmentID();
					 return (currentAssignmentID == appliedAssignment) && currentTasks.empty();
				}

				void setDispatcherStop()
				{
					m_dispatcher_running_flag.clear();
				}

				void setPersistentBufferStop()
				{
					m_persistentBuffer->StopPDB();
				}

				void setSkippedZero()
				{
					m_skipped.store(0);
				}
				
				std::atomic_int m_skipped;

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

				void ProcessAsyncAppend(const SizeType headID, int appendNum, std::string* appendPosting, std::function<void()> p_callback)
				{
					Append(headID, appendNum, appendPosting);

					if (p_callback != nullptr) {
						p_callback();
					}
					//clear running bit
					running[headID] = 0;
				}

				void ProcessAsyncDelete(const SizeType& p_id, std::function<void()> p_callback)
				{
					Delete(p_id);

					if (p_callback != nullptr) {
						p_callback();
					}
				}

				void ProcessAsyncPush(const Task& t)
				{
					// if (t.id == 0) {
					// 	LOG(Helper::LogLevel::LL_Info, "BUG: t.id = 0\n");
					// }
					while (!currentTasks.push(t)) {}
				}

				std::shared_ptr<VectorIndex> m_index;

				//std::unique_ptr<long long[]> m_vectorTranslateMap;
				//std::vector<long long> m_vectorTranslateMap;

				std::unique_ptr<IExtraSearcher<ValueType>> m_extraSearcher;

				//persisten buffer

				std::unique_ptr<PersistentBuffer> m_persistentBuffer;
				//std::unique_ptr<Dispatcher> m_dispatcher;

				std::unique_ptr<Helper::ThreadPool> m_threadPool;
				std::unique_ptr<Helper::ThreadPool> m_appendThreadPool;
				std::unique_ptr<Helper::ThreadPool> m_pushThreadPool;
				std::unique_ptr<Helper::ThreadPool> m_deleteThreadPool;
				std::unique_ptr<Helper::ThreadPool> m_pushPool;

				std::atomic<std::int32_t> m_tids;

				//for real-time service and release
				COMMON::Labelset m_deletedID;

				//delete mutex lock
				std::shared_timed_mutex m_dataDeleteLock;

				std::mutex m_dataAddLock;

				//config variable, set from beginning, will never change

				int m_internalResultNum;

				int m_resultNum;

				int m_replicaCount;

				int m_postingPageLimit;

				int m_postingVectorLimit;

				int m_vectorSize;

				int m_k = 2;

				bool m_clearHead = true;

				//atomic variable: m_vectornum for giving inserted vector a new id;

				std::atomic_int m_vectornum;

				//for search status
				
				int m_vectornum_last;

				int m_searchVectorLimit;

				int m_postingSize_avg;

				//insert information
				int m_split_num;

				//can't be atomic when append new entry, so need to set lock

				std::mutex m_headAddLock;
				
				//std::vector<std::unique_ptr<std::atomic_int>> m_postingSizes;
				std::vector<int> m_postingSizes;
				std::vector<float> m_postingRadius;

				int m_posting_num;

				//actually useless in AddHeadToPost = true
				SPTAG::COMMON::Dataset<long long> m_vectorTranslateMap;

				boost::lockfree::stack<ExtraWorkSpace*> m_workspaces;

				//dispatcher
				int appliedAssignment;
				int finishedAssignment;
				int m_appendThreadNum = 32;
				int m_deleteThreadNum = 1;
				int m_pushThreadNum = 5;
				std::atomic_flag m_dispatcher_running_flag;
				
				boost::lockfree::queue<Task, boost::lockfree::fixed_sized<true>> currentTasks;
				uint8_t* running;
				std::pair<SizeType, SizeType>* splitRoute;
				// std::unordered_map<SizeType, std::pair<SizeType, SizeType>> splitRoute;
			};
		}
	}
}
