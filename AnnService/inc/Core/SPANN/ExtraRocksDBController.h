// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#ifndef _SPTAG_SPANN_EXTRADBSEARCHER_H_
#define _SPTAG_SPANN_EXTRADBSEARCHER_H_

#include "inc/Helper/VectorSetReader.h"
#include "inc/Helper/AsyncFileReader.h"
#include "IExtraSearcher.h"
#include "ExtraFullGraphSearcher.h"
#include "../Common/TruthSet.h"
#include "inc/Helper/KeyValueIO.h"

#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "rocksdb/merge_operator.h"

#include <map>
#include <cmath>
#include <climits>
#include <future>

template <typename... Targs>
void DUMMY_CODE(Targs &&... /* unused */) {}

namespace SPTAG
{
    namespace SPANN
    {
        class RocksDBIO : public Helper::KeyValueIO
        {
        public:
            RocksDBIO() {}

            virtual ~RocksDBIO() { ShutDown(); }

            virtual bool Initialize(const char* filePath)
            {
                dbPath = std::move(std::string(filePath));
                dbOptions.create_if_missing = true;
                dbOptions.IncreaseParallelism();
                dbOptions.OptimizeLevelStyleCompaction();
                dbOptions.merge_operator.reset(new AnnMergeOperator);

                auto s = rocksdb::DB::Open(dbOptions, dbPath, &db);
                LOG(Helper::LogLevel::LL_Info, "SPFresh: New Rocksdb: %s\n", filePath);
                return s == rocksdb::Status::OK();
            }

            virtual void ShutDown() {
                db->Close();
                DestroyDB(dbPath, dbOptions);
                delete db;
            }

            virtual ErrorCode Get(const std::string& key, std::string* value) {
                auto s = db->Get(rocksdb::ReadOptions(), key, value);
                if (s == rocksdb::Status::OK()) {
                    return ErrorCode::Success;
                } else {
                    return ErrorCode::Fail;
                }
            }

            virtual ErrorCode Get(SizeType key, std::string* value) {
                return Get(Helper::Convert::Serialize<SizeType>(&key), value);
            }

            virtual ErrorCode Put(const std::string& key, const std::string& value) {
                auto s = db->Put(rocksdb::WriteOptions(), key, value);
                if (s == rocksdb::Status::OK()) {
                    return ErrorCode::Success;
                } else {
                    return ErrorCode::Fail;
                }
            }
            
            virtual ErrorCode Put(SizeType key, const std::string& value) {
                return Put(Helper::Convert::Serialize<SizeType>(&key), value);
            }

            virtual ErrorCode Put(SizeType key, SizeType id, const void* vector, SizeType dim) {
                using Helper::Convert::Serialize;
                std::string posting(std::move(Serialize<SizeType>(&id) + Serialize<SizeType>(vector, dim)));
                return Put(key, posting);
            }

            class AnnMergeOperator : public rocksdb::AssociativeMergeOperator
            {
            public:
                virtual bool Merge(const rocksdb::Slice& key, const rocksdb::Slice* existing_value,
                                   const rocksdb::Slice& value, std::string* new_value,
                                   rocksdb::Logger* logger) const override {
                    std::string newPosting;
                    if(existing_value)
                    {
                        newPosting += (*existing_value).ToString();
                        newPosting += value.ToString();
                    } else
                    {
                        newPosting += value.ToString();
                    }
                    *new_value = newPosting;
                    return true;
                }
                virtual const char* Name() const override {
                    return "AnnMergeOperator";
                }
            };

            ErrorCode Merge(SizeType key, const std::string& value) {
                if (value.size() == 0) {
                    LOG(Helper::LogLevel::LL_Error, "Error! empty append posting!\n");
                }
                auto s = db->Merge(rocksdb::WriteOptions(), Helper::Convert::Serialize<int>(&key, 1), value);
                if (s == rocksdb::Status::OK()) {
                    return ErrorCode::Success;
                } else {
                    return ErrorCode::Fail;
                }
            }

            ErrorCode Delete(SizeType key) {
                auto s = db->Delete(rocksdb::WriteOptions(), Helper::Convert::Serialize<int>(&key, 1));
                if (s == rocksdb::Status::OK()) {
                    return ErrorCode::Success;
                } else {
                    return ErrorCode::Fail;
                }
            }

            void ForceCompaction() {
                LOG(Helper::LogLevel::LL_Info, "Start Compaction\n");
                db->CompactRange(rocksdb::CompactRangeOptions(), nullptr, nullptr);
                LOG(Helper::LogLevel::LL_Info, "Finish Compaction\n");
            }
        private:
            std::string dbPath;
            rocksdb::DB* db;
            rocksdb::Options dbOptions;
        };
 
        template <typename ValueType>
        class ExtraRocksDBController : public IExtraSearcher
        {
        private:
            RocksDBIO db;
            std::atomic_uint64_t m_postingNum;
        public:
            ExtraRocksDBController(const char* dbPath, int dim) { db.Initialize(dbPath); m_vectorInfoSize = dim * sizeof(ValueType) + sizeof(int);}

            virtual ~ExtraRocksDBController() { }

            virtual bool LoadIndex(Options& p_opt) {
                /*
                m_extraFullGraphFile = p_opt.m_indexDirectory + FolderSep + p_opt.m_ssdIndex;
                std::string curFile = m_extraFullGraphFile;
                do {
                    auto curIndexFile = f_createAsyncIO();
                    if (curIndexFile == nullptr || !curIndexFile->Initialize(curFile.c_str(), std::ios::binary | std::ios::in, 
#ifdef BATCH_READ
                        p_opt.m_searchInternalResultNum, 2, 2, p_opt.m_iSSDNumberOfThreads
#else
                        p_opt.m_searchInternalResultNum * p_opt.m_iSSDNumberOfThreads / p_opt.m_ioThreads + 1, 2, 2, p_opt.m_ioThreads
#endif
                    )) {
                        LOG(Helper::LogLevel::LL_Error, "Cannot open file:%s!\n", curFile.c_str());
                        return false;
                    }

                    m_indexFiles.emplace_back(curIndexFile);
                    m_listInfos.emplace_back(0);
                    m_totalListCount += LoadingHeadInfo(curFile, p_opt.m_searchPostingPageLimit, m_listInfos.back());

                    curFile = m_extraFullGraphFile + "_" + std::to_string(m_indexFiles.size());
                } while (fileexists(curFile.c_str()));
                m_listPerFile = static_cast<int>((m_totalListCount + m_indexFiles.size() - 1) / m_indexFiles.size());

#ifndef _MSC_VER
                Helper::AIOTimeout.tv_nsec = p_opt.m_iotimeout * 1000;
#endif
                */
                return true;
            }

            virtual void SearchIndex(ExtraWorkSpace* p_exWorkSpace,
                QueryResult& p_queryResults,
                std::shared_ptr<VectorIndex> p_index,
                SearchStats* p_stats, std::set<int>* truth, std::map<int, std::set<int>>* found)
            {
                const auto postingListCount = static_cast<uint32_t>(p_exWorkSpace->m_postingIDs.size());

                p_exWorkSpace->m_deduper.clear();

                COMMON::QueryResultSet<ValueType>& queryResults = *((COMMON::QueryResultSet<ValueType>*)&p_queryResults);

                int diskRead = 0;
                int diskIO = 0;
                int listElements = 0;

                for (uint32_t pi = 0; pi < postingListCount; ++pi)
                {
                    auto curPostingID = p_exWorkSpace->m_postingIDs[pi];
                    auto postingP = new std::string;

                    SearchIndex(curPostingID, *postingP);
                    p_exWorkSpace->m_pageBuffers[pi].SetPointer(std::shared_ptr<uint8_t>(
                            (uint8_t *)&(*postingP)[0], [postingP](uint8_t*){ delete postingP; }));
                    int vectorNum = postingP->size() / m_vectorInfoSize;

                    diskIO++;
                    diskRead++;
                    listElements += vectorNum;
 
                    auto buffer = reinterpret_cast<char*>((p_exWorkSpace->m_pageBuffers[pi]).GetBuffer());
                    for (int i = 0; i < vectorNum; i++) {
                        char* vectorInfo = postingP->data() + i * m_vectorInfoSize;
                        int vectorID = *(reinterpret_cast<int*>(vectorInfo));
                        if (p_exWorkSpace->m_deduper.CheckAndSet(vectorID)) continue;
                        auto distance2leaf = p_index->ComputeDistance(queryResults.GetQuantizedTarget(), vectorInfo + sizeof(int));
                        queryResults.AddPoint(vectorID, distance2leaf);
                    }
                    if (truth) {
                        for (int i = 0; i < vectorNum; ++i) {
                            char* vectorInfo = buffer + i * m_vectorInfoSize;
                            int vectorID = *(reinterpret_cast<int*>(vectorInfo));
                            if (truth && truth->count(vectorID))
                                (*found)[curPostingID].insert(vectorID);
                        }
                    }
                }

                if (p_stats)
                {
                    p_stats->m_totalListElementsCount = listElements;
                    p_stats->m_diskIOCount = diskIO;
                    p_stats->m_diskAccessCount = diskRead;
                }
            }

            bool BuildIndex(std::shared_ptr<Helper::VectorSetReader>& p_reader, std::shared_ptr<VectorIndex> p_headIndex, Options& p_opt) {
                std::string outputFile = p_opt.m_indexDirectory + FolderSep + p_opt.m_ssdIndex;
                if (outputFile.empty())
                {
                    LOG(Helper::LogLevel::LL_Error, "Output file can't be empty!\n");
                    return false;
                }

                int numThreads = p_opt.m_iSSDNumberOfThreads;
                int candidateNum = p_opt.m_internalResultNum;

                std::unordered_set<SizeType> headVectorIDS;
                if (p_opt.m_headIDFile.empty()) {
                    LOG(Helper::LogLevel::LL_Error, "Not found VectorIDTranslate!\n");
                    return false;
                }

                {
                    auto ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize((p_opt.m_indexDirectory + FolderSep +  p_opt.m_headIDFile).c_str(), std::ios::binary | std::ios::in)) {
                        LOG(Helper::LogLevel::LL_Error, "failed open VectorIDTranslate: %s\n", p_opt.m_headIDFile.c_str());
                        return false;
                    }

                    std::uint64_t vid;
                    while (ptr->ReadBinary(sizeof(vid), reinterpret_cast<char*>(&vid)) == sizeof(vid))
                    {
                        headVectorIDS.insert(static_cast<SizeType>(vid));
                    }
                    LOG(Helper::LogLevel::LL_Info, "Loaded %u Vector IDs\n", static_cast<uint32_t>(headVectorIDS.size()));
                }

                SizeType fullCount = 0;
                size_t vectorInfoSize = 0;
                {
                    auto fullVectors = p_reader->GetVectorSet();
                    fullCount = fullVectors->Count();
                    vectorInfoSize = fullVectors->PerVectorDataSize() + sizeof(int);
                }

                Selection selections(static_cast<size_t>(fullCount) * p_opt.m_replicaCount, p_opt.m_tmpdir);
                LOG(Helper::LogLevel::LL_Info, "Full vector count:%d Edge bytes:%llu selection size:%zu, capacity size:%zu\n", fullCount, sizeof(Edge), selections.m_selections.size(), selections.m_selections.capacity());
                std::vector<std::atomic_int> replicaCount(fullCount);
                std::vector<std::atomic_int> postingListSize(headVectorIDS.size());
                for (auto& pls : postingListSize) pls = 0;
                std::unordered_set<SizeType> emptySet;
                SizeType batchSize = (fullCount + p_opt.m_batches - 1) / p_opt.m_batches;

                auto t1 = std::chrono::high_resolution_clock::now();
                if (p_opt.m_batches > 1) selections.SaveBatch();
                {
                    LOG(Helper::LogLevel::LL_Info, "Preparation done, start candidate searching.\n");
                    SizeType sampleSize = p_opt.m_samples;
                    std::vector<SizeType> samples(sampleSize, 0);
                    for (int i = 0; i < p_opt.m_batches; i++) {
                        SizeType start = i * batchSize;
                        SizeType end = min(start + batchSize, fullCount);
                        auto fullVectors = p_reader->GetVectorSet(start, end);
                        if (p_opt.m_distCalcMethod == DistCalcMethod::Cosine && !p_reader->IsNormalized()) fullVectors->Normalize(p_opt.m_iSSDNumberOfThreads);

                        emptySet.clear();

                        int sampleNum = 0;
                        for (int j = start; j < end && sampleNum < sampleSize; j++)
                        {
                            if (headVectorIDS.count(j) == 0) samples[sampleNum++] = j - start;
                        }

                        float acc = 0;
#pragma omp parallel for schedule(dynamic)
                        for (int j = 0; j < sampleNum; j++)
                        {
                            COMMON::Utils::atomic_float_add(&acc, COMMON::TruthSet::CalculateRecall(p_headIndex.get(), fullVectors->GetVector(samples[j]), candidateNum));
                        }
                        acc = acc / sampleNum;
                        LOG(Helper::LogLevel::LL_Info, "Batch %d vector(%d,%d) loaded with %d vectors (%zu) HeadIndex acc @%d:%f.\n", i, start, end, fullVectors->Count(), selections.m_selections.size(), candidateNum, acc);

                        p_headIndex->ApproximateRNG(fullVectors, emptySet, candidateNum, selections.m_selections.data(), p_opt.m_replicaCount, numThreads, p_opt.m_gpuSSDNumTrees, p_opt.m_gpuSSDLeafSize, p_opt.m_rngFactor, p_opt.m_numGPUs);

                        for (SizeType j = start; j < end; j++) {
                            replicaCount[j] = 0;
                            size_t vecOffset = j * (size_t)p_opt.m_replicaCount;
                            for (int resNum = 0; resNum < p_opt.m_replicaCount && selections[vecOffset + resNum].node != INT_MAX; resNum++) {
                                ++postingListSize[selections[vecOffset + resNum].node];
                                selections[vecOffset + resNum].tonode = j;
                                ++replicaCount[j];
                            }
                        }

                        if (p_opt.m_batches > 1) selections.SaveBatch();
                    }
                }
                auto t2 = std::chrono::high_resolution_clock::now();
                LOG(Helper::LogLevel::LL_Info, "Searching replicas ended. Search Time: %.2lf mins\n", ((double)std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count()) / 60.0);

                if (p_opt.m_batches > 1) selections.LoadBatch(0, static_cast<size_t>(fullCount) * p_opt.m_replicaCount);

                // Sort results either in CPU or GPU
                VectorIndex::SortSelections(&selections.m_selections);

                auto t3 = std::chrono::high_resolution_clock::now();
                LOG(Helper::LogLevel::LL_Info, "Time to sort selections:%.2lf sec.\n", ((double)std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count()) + ((double)std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count()) / 1000);

                int postingSizeLimit = INT_MAX;
                if (p_opt.m_postingPageLimit > 0)
                {
                    postingSizeLimit = static_cast<int>(p_opt.m_postingPageLimit * PageSize / vectorInfoSize);
                }

                LOG(Helper::LogLevel::LL_Info, "Posting size limit: %d\n", postingSizeLimit);

                {
                    std::vector<int> replicaCountDist(p_opt.m_replicaCount + 1, 0);
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {
                        ++replicaCountDist[replicaCount[i]];
                    }

                    LOG(Helper::LogLevel::LL_Info, "Before Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        LOG(Helper::LogLevel::LL_Info, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
                    }
                }

#pragma omp parallel for schedule(dynamic)
                for (int i = 0; i < postingListSize.size(); ++i)
                {
                    if (postingListSize[i] <= postingSizeLimit) continue;

                    std::size_t selectIdx = std::lower_bound(selections.m_selections.begin(), selections.m_selections.end(), i, Selection::g_edgeComparer) - selections.m_selections.begin();

                    for (size_t dropID = postingSizeLimit; dropID < postingListSize[i]; ++dropID)
                    {
                        int tonode = selections.m_selections[selectIdx + dropID].tonode;
                        --replicaCount[tonode];
                    }
                    postingListSize[i] = postingSizeLimit;
                }

                if (p_opt.m_outputEmptyReplicaID)
                {
                    std::vector<int> replicaCountDist(p_opt.m_replicaCount + 1, 0);
                    auto ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize("EmptyReplicaID.bin", std::ios::binary | std::ios::out)) {
                        LOG(Helper::LogLevel::LL_Error, "Fail to create EmptyReplicaID.bin!\n");
                        return false;
                    }
                    for (int i = 0; i < replicaCount.size(); ++i)
                    {

                        ++replicaCountDist[replicaCount[i]];

                        if (replicaCount[i] < 2)
                        {
                            long long vid = i;
                            if (ptr->WriteBinary(sizeof(vid), reinterpret_cast<char*>(&vid)) != sizeof(vid)) {
                                LOG(Helper::LogLevel::LL_Error, "Failt to write EmptyReplicaID.bin!");
                                return false;
                            }
                        }
                    }

                    LOG(Helper::LogLevel::LL_Info, "After Posting Cut:\n");
                    for (int i = 0; i < replicaCountDist.size(); ++i)
                    {
                        LOG(Helper::LogLevel::LL_Info, "Replica Count Dist: %d, %d\n", i, replicaCountDist[i]);
                    }
                }

                auto t4 = std::chrono::high_resolution_clock::now();
                LOG(SPTAG::Helper::LogLevel::LL_Info, "Time to perform posting cut:%.2lf sec.\n", ((double)std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count()) + ((double)std::chrono::duration_cast<std::chrono::milliseconds>(t4 - t3).count()) / 1000);

                if (p_opt.m_ssdIndexFileNum > 1) selections.SaveBatch();

                auto fullVectors = p_reader->GetVectorSet();
                if (p_opt.m_distCalcMethod == DistCalcMethod::Cosine && !p_reader->IsNormalized()) fullVectors->Normalize(p_opt.m_iSSDNumberOfThreads);

                for (int id = 0; id < postingListSize.size(); id++) 
                {
                    std::string postinglist;
                    std::size_t selectIdx = selections.lower_bound(id);
                    for (int j = 0; j < postingListSize[id]; ++j) {
                        if (selections[selectIdx].node != id) {
                            LOG(Helper::LogLevel::LL_Error, "Selection ID NOT MATCH\n");
                            exit(1);
                        }
                        int fullID = selections[selectIdx++].tonode;
                        size_t dim = fullVectors->Dimension();
                        // First Vector ID, then Vector
                        postinglist += Helper::Convert::Serialize<int>(&fullID, 1);
                        postinglist += Helper::Convert::Serialize<ValueType>(fullVectors->GetVector(fullID), dim);
                    }

                    AddIndex(id, postinglist);
                }

                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(p_opt.m_ssdInfoFile.c_str(), std::ios::binary | std::ios::out)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed open file %s\n", p_opt.m_ssdInfoFile.c_str());
                    exit(1);
                }
                //Number of all documents.
                int i32Val = static_cast<int>(fullCount);
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
                    i32Val = postingListSize[id].load();
                    if (ptr->WriteBinary(sizeof(i32Val), reinterpret_cast<char*>(&i32Val)) != sizeof(i32Val)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to write SSDIndexInfo File!");
                        exit(1);
                    }
                }

                LOG(Helper::LogLevel::LL_Info, "SPFresh: initialize deleteMap\n");
                COMMON::Labelset m_deleteID;
                m_deleteID.Initialize(fullCount, p_headIndex->m_iDataBlockSize, p_headIndex->m_iDataCapacity);
                LOG(Helper::LogLevel::LL_Info, "SPFresh: save deleteMap\n");
                m_deleteID.Save(p_opt.m_fullDeletedIDFile);

                auto t5 = std::chrono::high_resolution_clock::now();
                double elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(t5 - t1).count();
                LOG(Helper::LogLevel::LL_Info, "Total used time: %.2lf minutes (about %.2lf hours).\n", elapsedSeconds / 60.0, elapsedSeconds / 3600.0);
            
                return true;
            }

            virtual ErrorCode AppendPosting(SizeType headID, const std::string& appendPosting) {
                if (appendPosting.size() == 0) {
                    LOG(Helper::LogLevel::LL_Error, "Error! empty append posting!\n");
                }
                return db.Merge(headID, appendPosting);
            }

            virtual void ForceCompaction() {db.ForceCompaction();}

            virtual ErrorCode SearchIndex(SizeType headID, std::string& posting) {  return db.Get(headID, &posting); }
            virtual ErrorCode AddIndex(SizeType headID, const std::string& posting) { m_postingNum++; return db.Put(headID, posting); }
            virtual ErrorCode DeleteIndex(SizeType headID) { m_postingNum--; return db.Delete(headID); }
            virtual SizeType  GetIndexSize() { return m_postingNum; }
        private:
            struct ListInfo
            {
                int listEleCount = 0;

                std::uint16_t listPageCount = 0;

                std::uint64_t listOffset = 0;

                std::uint16_t pageOffset = 0;
            };

            int LoadingHeadInfo(const std::string& p_file, int p_postingPageLimit, std::vector<ListInfo>& m_listInfos)
            {
                auto ptr = SPTAG::f_createIO();
                if (ptr == nullptr || !ptr->Initialize(p_file.c_str(), std::ios::binary | std::ios::in)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to open file: %s\n", p_file.c_str());
                    exit(1);
                }

                int m_listCount;
                int m_totalDocumentCount;
                int m_iDataDimension;
                int m_listPageOffset;

                if (ptr->ReadBinary(sizeof(m_listCount), reinterpret_cast<char*>(&m_listCount)) != sizeof(m_listCount)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                    exit(1);
                }
                if (ptr->ReadBinary(sizeof(m_totalDocumentCount), reinterpret_cast<char*>(&m_totalDocumentCount)) != sizeof(m_totalDocumentCount)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                    exit(1);
                }
                if (ptr->ReadBinary(sizeof(m_iDataDimension), reinterpret_cast<char*>(&m_iDataDimension)) != sizeof(m_iDataDimension)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                    exit(1);
                }
                if (ptr->ReadBinary(sizeof(m_listPageOffset), reinterpret_cast<char*>(&m_listPageOffset)) != sizeof(m_listPageOffset)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                    exit(1);
                }

                if (m_vectorInfoSize == 0) m_vectorInfoSize = m_iDataDimension * sizeof(ValueType) + sizeof(int);
                else if (m_vectorInfoSize != m_iDataDimension * sizeof(ValueType) + sizeof(int)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head info file! DataDimension and ValueType are not match!\n");
                    exit(1);
                }

                m_listInfos.resize(m_listCount);

                size_t totalListElementCount = 0;

                std::map<int, int> pageCountDist;

                size_t biglistCount = 0;
                size_t biglistElementCount = 0;
                int pageNum;
                for (int i = 0; i < m_listCount; ++i)
                {
                    if (ptr->ReadBinary(sizeof(pageNum), reinterpret_cast<char*>(&(pageNum))) != sizeof(pageNum)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                        exit(1);
                    }
                    if (ptr->ReadBinary(sizeof(m_listInfos[i].pageOffset), reinterpret_cast<char*>(&(m_listInfos[i].pageOffset))) != sizeof(m_listInfos[i].pageOffset)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                        exit(1);
                    }
                    if (ptr->ReadBinary(sizeof(m_listInfos[i].listEleCount), reinterpret_cast<char*>(&(m_listInfos[i].listEleCount))) != sizeof(m_listInfos[i].listEleCount)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                        exit(1);
                    }
                    if (ptr->ReadBinary(sizeof(m_listInfos[i].listPageCount), reinterpret_cast<char*>(&(m_listInfos[i].listPageCount))) != sizeof(m_listInfos[i].listPageCount)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to read head info file!\n");
                        exit(1);
                    }

                    m_listInfos[i].listOffset = (static_cast<uint64_t>(m_listPageOffset + pageNum) << PageSizeEx);
                    m_listInfos[i].listEleCount = min(m_listInfos[i].listEleCount, (min(static_cast<int>(m_listInfos[i].listPageCount), p_postingPageLimit) << PageSizeEx) / m_vectorInfoSize);
                    m_listInfos[i].listPageCount = static_cast<std::uint16_t>(ceil((m_vectorInfoSize * m_listInfos[i].listEleCount + m_listInfos[i].pageOffset) * 1.0 / (1 << PageSizeEx)));
                    totalListElementCount += m_listInfos[i].listEleCount;
                    int pageCount = m_listInfos[i].listPageCount;

                    if (pageCount > 1)
                    {
                        ++biglistCount;
                        biglistElementCount += m_listInfos[i].listEleCount;
                    }

                    if (pageCountDist.count(pageCount) == 0)
                    {
                        pageCountDist[pageCount] = 1;
                    }
                    else
                    {
                        pageCountDist[pageCount] += 1;
                    }
                }

                LOG(Helper::LogLevel::LL_Info,
                    "Finish reading header info, list count %d, total doc count %d, dimension %d, list page offset %d.\n",
                    m_listCount,
                    m_totalDocumentCount,
                    m_iDataDimension,
                    m_listPageOffset);


                LOG(Helper::LogLevel::LL_Info,
                    "Big page (>4K): list count %zu, total element count %zu.\n",
                    biglistCount,
                    biglistElementCount);

                LOG(Helper::LogLevel::LL_Info, "Total Element Count: %llu\n", totalListElementCount);

                for (auto& ele : pageCountDist)
                {
                    LOG(Helper::LogLevel::LL_Info, "Page Count Dist: %d %d\n", ele.first, ele.second);
                }

                return m_listCount;
            }

        private:
            
            std::string m_extraFullGraphFile;

            std::vector<std::vector<ListInfo>> m_listInfos;

            std::vector<std::shared_ptr<Helper::DiskPriorityIO>> m_indexFiles;

            int m_vectorInfoSize = 0;

            int m_totalListCount = 0;

            int m_listPerFile = 0;
        };
    } // namespace SPANN
} // namespace SPTAG

#endif // _SPTAG_SPANN_EXTRADBSEARCHER_H_
