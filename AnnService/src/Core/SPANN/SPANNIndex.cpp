// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#include "inc/Core/SPANN/Index.h"
#include "inc/Helper/VectorSetReaders/MemoryReader.h"
#include "inc/Core/SPANN/ExtraFullGraphSearcher.h"
#include "inc/Core/SPANN/ExtraRocksDBController.h"
#include <shared_mutex>
#include <chrono>
#include <random>

#pragma warning(disable:4242)  // '=' : conversion from 'int' to 'short', possible loss of data
#pragma warning(disable:4244)  // '=' : conversion from 'int' to 'short', possible loss of data
#pragma warning(disable:4127)  // conditional expression is constant

namespace SPTAG
{
    namespace SPANN
    {
        std::atomic_int ExtraWorkSpace::g_spaceCount(0);
        EdgeCompare Selection::g_edgeComparer;

        std::function<std::shared_ptr<Helper::DiskPriorityIO>(void)> f_createAsyncIO = []() -> std::shared_ptr<Helper::DiskPriorityIO> { return std::shared_ptr<Helper::DiskPriorityIO>(new Helper::AsyncFileIO()); };

        template <typename T>
        bool Index<T>::CheckHeadIndexType() {
            SPTAG::VectorValueType v1 = m_index->GetVectorValueType(), v2 = GetEnumValueType<T>();
            if (v1 != v2) {
                LOG(Helper::LogLevel::LL_Error, "Head index and vectors don't have the same value types, which are %s %s\n",
                    SPTAG::Helper::Convert::ConvertToString(v1).c_str(),
                    SPTAG::Helper::Convert::ConvertToString(v2).c_str()
                );
                if (!SPTAG::COMMON::DistanceUtils::Quantizer) return false;
            }
            return true;
        }

        template <typename T>
        ErrorCode Index<T>::LoadConfig(Helper::IniReader& p_reader)
        {
            IndexAlgoType algoType = p_reader.GetParameter("Base", "IndexAlgoType", IndexAlgoType::Undefined);
            VectorValueType valueType = p_reader.GetParameter("Base", "ValueType", VectorValueType::Undefined);
            if ((m_index = CreateInstance(algoType, valueType)) == nullptr) return ErrorCode::FailedParseValue;

            std::string sections[] = { "Base", "SelectHead", "BuildHead", "BuildSSDIndex" };
            for (int i = 0; i < 4; i++) {
                auto parameters = p_reader.GetParameters(sections[i].c_str());
                for (auto iter = parameters.begin(); iter != parameters.end(); iter++) {
                    SetParameter(iter->first.c_str(), iter->second.c_str(), sections[i].c_str());
                }
            }
            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::LoadIndexDataFromMemory(const std::vector<ByteArray>& p_indexBlobs)
        {
            if (m_index->LoadIndexDataFromMemory(p_indexBlobs) != ErrorCode::Success) return ErrorCode::Fail;

            m_index->SetParameter("NumberOfThreads", std::to_string(m_options.m_iSSDNumberOfThreads));
            m_index->SetParameter("MaxCheck", std::to_string(m_options.m_maxCheck));
            m_index->SetParameter("HashTableExponent", std::to_string(m_options.m_hashExp));
            m_index->UpdateIndex();
            m_index->SetReady(true);

            m_extraSearcher.reset(new ExtraFullGraphSearcher<T>());
            if (!m_extraSearcher->LoadIndex(m_options)) return ErrorCode::Fail;

            m_vectorTranslateMap.reset((std::uint64_t*)(p_indexBlobs.back().Data()), [=](std::uint64_t* ptr) {});

            omp_set_num_threads(m_options.m_iSSDNumberOfThreads);
            m_workSpacePool.reset(new COMMON::WorkSpacePool<ExtraWorkSpace>());
            m_workSpacePool->Init(m_options.m_iSSDNumberOfThreads, m_options.m_maxCheck, m_options.m_hashExp, m_options.m_searchInternalResultNum, min(m_options.m_postingPageLimit, m_options.m_searchPostingPageLimit + 1) << PageSizeEx);
            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::LoadIndexData(const std::vector<std::shared_ptr<Helper::DiskPriorityIO>>& p_indexStreams)
        {
            if (m_index->LoadIndexData(p_indexStreams) != ErrorCode::Success) return ErrorCode::Fail;

            m_index->SetParameter("NumberOfThreads", std::to_string(m_options.m_iSSDNumberOfThreads));
            m_index->SetParameter("MaxCheck", std::to_string(m_options.m_maxCheck));
            m_index->SetParameter("HashTableExponent", std::to_string(m_options.m_hashExp));
            m_index->UpdateIndex();
            m_index->SetReady(true);

            // TODO: Choose an extra searcher based on config
            m_extraSearcher.reset(new ExtraFullGraphSearcher<T>());
            if (!m_extraSearcher->LoadIndex(m_options)) return ErrorCode::Fail;

            m_vectorTranslateMap.reset(new std::uint64_t[m_index->GetNumSamples()], std::default_delete<std::uint64_t[]>());
            IOBINARY(p_indexStreams.back(), ReadBinary, sizeof(std::uint64_t) * m_index->GetNumSamples(), reinterpret_cast<char*>(m_vectorTranslateMap.get()));

            omp_set_num_threads(m_options.m_iSSDNumberOfThreads);
            m_workSpacePool = std::make_unique<COMMON::WorkSpacePool<ExtraWorkSpace>>();
            m_workSpacePool->Init(m_options.m_iSSDNumberOfThreads, m_options.m_maxCheck, m_options.m_hashExp, m_options.m_searchInternalResultNum, min(m_options.m_postingPageLimit, m_options.m_searchPostingPageLimit + 1) << PageSizeEx);

            m_deletedID.Load(m_options.m_fullDeletedIDFile, m_iDataBlockSize, m_iDataCapacity);

            // TODO: choose a proper size
            m_rwLocks = std::make_unique<std::shared_timed_mutex[]>(500000000);
            m_postingSizes = std::make_unique<std::atomic_uint32_t[]>(500000000);

            for (int idx = 0; idx < m_extraSearcher->GetIndexSize(); idx++) {
                uint32_t tmp;
                IOBINARY(p_indexStreams.back(), ReadBinary, sizeof(uint32_t), reinterpret_cast<char*>(&tmp));
                m_postingSizes[idx].store(tmp);
            }

            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::SaveConfig(std::shared_ptr<Helper::DiskPriorityIO> p_configOut)
        {
            IOSTRING(p_configOut, WriteString, "[Base]\n");
#define DefineBasicParameter(VarName, VarType, DefaultValue, RepresentStr) \
                IOSTRING(p_configOut, WriteString, (RepresentStr + std::string("=") + SPTAG::Helper::Convert::ConvertToString(m_options.VarName) + std::string("\n")).c_str()); \

#include "inc/Core/SPANN/ParameterDefinitionList.h"
#undef DefineBasicParameter

            IOSTRING(p_configOut, WriteString, "[SelectHead]\n");
#define DefineSelectHeadParameter(VarName, VarType, DefaultValue, RepresentStr) \
                IOSTRING(p_configOut, WriteString, (RepresentStr + std::string("=") + SPTAG::Helper::Convert::ConvertToString(m_options.VarName) + std::string("\n")).c_str()); \

#include "inc/Core/SPANN/ParameterDefinitionList.h"
#undef DefineSelectHeadParameter

            IOSTRING(p_configOut, WriteString, "[BuildHead]\n");
#define DefineBuildHeadParameter(VarName, VarType, DefaultValue, RepresentStr) \
                IOSTRING(p_configOut, WriteString, (RepresentStr + std::string("=") + SPTAG::Helper::Convert::ConvertToString(m_options.VarName) + std::string("\n")).c_str()); \

#include "inc/Core/SPANN/ParameterDefinitionList.h"
#undef DefineBuildHeadParameter

            m_index->SaveConfig(p_configOut);

            Helper::Convert::ConvertStringTo<int>(m_index->GetParameter("HashTableExponent").c_str(), m_options.m_hashExp);
            IOSTRING(p_configOut, WriteString, "[BuildSSDIndex]\n");
#define DefineSSDParameter(VarName, VarType, DefaultValue, RepresentStr) \
                IOSTRING(p_configOut, WriteString, (RepresentStr + std::string("=") + SPTAG::Helper::Convert::ConvertToString(m_options.VarName) + std::string("\n")).c_str()); \

#include "inc/Core/SPANN/ParameterDefinitionList.h"
#undef DefineSSDParameter

            IOSTRING(p_configOut, WriteString, "\n");
            return ErrorCode::Success;
        }

        template<typename T>
        ErrorCode Index<T>::SaveIndexData(const std::vector<std::shared_ptr<Helper::DiskPriorityIO>>& p_indexStreams)
        {
            if (m_index == nullptr || m_vectorTranslateMap == nullptr) return ErrorCode::EmptyIndex;

            ErrorCode ret;
            if ((ret = m_index->SaveIndexData(p_indexStreams)) != ErrorCode::Success) return ret;

            IOBINARY(p_indexStreams.back(), WriteBinary, sizeof(std::uint64_t) * m_index->GetNumSamples(), (char*)(m_vectorTranslateMap.get()));
            m_deletedID.Save(m_options.m_fullDeletedIDFile);
            return ErrorCode::Success;
        }

#pragma region K-NN search

        template<typename T>
        ErrorCode Index<T>::SearchIndex(QueryResult &p_query, bool p_searchDeleted) const
        {
            if (!m_bReady) return ErrorCode::EmptyIndex;

            m_index->SearchIndex(p_query);

            auto* p_queryResults = (COMMON::QueryResultSet<T>*) & p_query;
            std::shared_ptr<ExtraWorkSpace> workSpace = nullptr;
            if (m_extraSearcher != nullptr) {
                workSpace = m_workSpacePool->Rent();
                workSpace->m_postingIDs.clear();

                float limitDist = p_queryResults->GetResult(0)->Dist * m_options.m_maxDistRatio;
                for (int i = 0; i < m_options.m_searchInternalResultNum; ++i)
                {
                    auto res = p_queryResults->GetResult(i);
                    if (res->VID == -1 || (limitDist > 0.1 && res->Dist > limitDist)) break;
                    workSpace->m_postingIDs.emplace_back(res->VID);
                }

                for (int i = 0; i < p_queryResults->GetResultNum(); ++i)
                {
                    auto res = p_queryResults->GetResult(i);
                    if (res->VID == -1) break;
                    res->VID = static_cast<SizeType>((m_vectorTranslateMap.get())[res->VID]);
                }

                p_queryResults->Reverse();
                m_extraSearcher->SearchIndex(workSpace.get(), *p_queryResults, m_index, nullptr);
                p_queryResults->SortResult();
                m_workSpacePool->Return(workSpace);
            }

            if (p_query.WithMeta() && nullptr != m_pMetadata)
            {
                for (int i = 0; i < p_query.GetResultNum(); ++i)
                {
                    SizeType result = p_query.GetResult(i)->VID;
                    p_query.SetMetadata(i, (result < 0) ? ByteArray::c_empty : m_pMetadata->GetMetadataCopy(result));
                }
            }
            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::DebugSearchDiskIndex(QueryResult& p_query, int p_subInternalResultNum, int p_internalResultNum,
                                                 SearchStats* p_stats, std::set<int>* truth, std::map<int, std::set<int>>* found)
        {
            if (nullptr == m_extraSearcher) return ErrorCode::EmptyIndex;

            COMMON::QueryResultSet<T> newResults(*((COMMON::QueryResultSet<T>*)&p_query));
            for (int i = 0; i < newResults.GetResultNum(); ++i)
            {
                auto res = newResults.GetResult(i);
                if (res->VID == -1) break;

                auto global_VID = static_cast<SizeType>((m_vectorTranslateMap.get())[res->VID]);
                if (truth && truth->count(global_VID)) (*found)[res->VID].insert(global_VID);
                res->VID = global_VID;
            }
            newResults.Reverse();

            auto auto_ws = m_workSpacePool->Rent();

            int partitions = (p_internalResultNum + p_subInternalResultNum - 1) / p_subInternalResultNum;
            float limitDist = p_query.GetResult(0)->Dist * m_options.m_maxDistRatio;
            for (SizeType p = 0; p < partitions; p++) {
                int subInternalResultNum = min(p_subInternalResultNum, p_internalResultNum - p_subInternalResultNum * p);

                auto_ws->m_postingIDs.clear();

                for (int i = p * p_subInternalResultNum; i < p * p_subInternalResultNum + subInternalResultNum; i++)
                {
                    auto res = p_query.GetResult(i);
                    if (res->VID == -1 || (limitDist > 0.1 && res->Dist > limitDist)) break;
                    auto_ws->m_postingIDs.emplace_back(res->VID);
                }

                m_extraSearcher->SearchIndex(auto_ws.get(), newResults, m_index, p_stats, truth, found);
            }

            m_workSpacePool->Return(auto_ws);

            newResults.SortResult();
            std::copy(newResults.GetResults(), newResults.GetResults() + newResults.GetResultNum(), p_query.GetResults());
            return ErrorCode::Success;
        }
#pragma endregion

        template <typename T>
        void Index<T>::SelectHeadAdjustOptions(int p_vectorCount) {
            if (m_options.m_headVectorCount != 0) m_options.m_ratio = m_options.m_headVectorCount * 1.0 / p_vectorCount;
            int headCnt = static_cast<int>(std::round(m_options.m_ratio * p_vectorCount));
            if (headCnt == 0)
            {
                for (double minCnt = 1; headCnt == 0; minCnt += 0.2)
                {
                    m_options.m_ratio = minCnt / p_vectorCount;
                    headCnt = static_cast<int>(std::round(m_options.m_ratio * p_vectorCount));
                }

                LOG(Helper::LogLevel::LL_Info, "Setting requires to select none vectors as head, adjusted it to %d vectors\n", headCnt);
            }

            if (m_options.m_iBKTKmeansK > headCnt)
            {
                m_options.m_iBKTKmeansK = headCnt;
                LOG(Helper::LogLevel::LL_Info, "Setting of cluster number is less than head count, adjust it to %d\n", headCnt);
            }

            if (m_options.m_selectThreshold == 0)
            {
                m_options.m_selectThreshold = min(p_vectorCount - 1, static_cast<int>(1 / m_options.m_ratio));
                LOG(Helper::LogLevel::LL_Info, "Set SelectThreshold to %d\n", m_options.m_selectThreshold);
            }

            if (m_options.m_splitThreshold == 0)
            {
                m_options.m_splitThreshold = min(p_vectorCount - 1, static_cast<int>(m_options.m_selectThreshold * 2));
                LOG(Helper::LogLevel::LL_Info, "Set SplitThreshold to %d\n", m_options.m_splitThreshold);
            }

            if (m_options.m_splitFactor == 0)
            {
                m_options.m_splitFactor = min(p_vectorCount - 1, static_cast<int>(std::round(1 / m_options.m_ratio) + 0.5));
                LOG(Helper::LogLevel::LL_Info, "Set SplitFactor to %d\n", m_options.m_splitFactor);
            }
        }

        template <typename T>
        int Index<T>::SelectHeadDynamicallyInternal(const std::shared_ptr<COMMON::BKTree> p_tree, int p_nodeID,
                                                    const Options& p_opts, std::vector<int>& p_selected)
        {
            typedef std::pair<int, int> CSPair;
            std::vector<CSPair> children;
            int childrenSize = 1;
            const auto& node = (*p_tree)[p_nodeID];
            if (node.childStart >= 0)
            {
                children.reserve(node.childEnd - node.childStart);
                for (int i = node.childStart; i < node.childEnd; ++i)
                {
                    int cs = SelectHeadDynamicallyInternal(p_tree, i, p_opts, p_selected);
                    if (cs > 0)
                    {
                        children.emplace_back(i, cs);
                        childrenSize += cs;
                    }
                }
            }

            if (childrenSize >= p_opts.m_selectThreshold)
            {
                if (node.centerid < (*p_tree)[0].centerid)
                {
                    p_selected.push_back(node.centerid);
                }

                if (childrenSize > p_opts.m_splitThreshold)
                {
                    std::sort(children.begin(), children.end(), [](const CSPair& a, const CSPair& b)
                    {
                        return a.second > b.second;
                    });

                    size_t selectCnt = static_cast<size_t>(std::ceil(childrenSize * 1.0 / p_opts.m_splitFactor) + 0.5);
                    //if (selectCnt > 1) selectCnt -= 1;
                    for (size_t i = 0; i < selectCnt && i < children.size(); ++i)
                    {
                        p_selected.push_back((*p_tree)[children[i].first].centerid);
                    }
                }

                return 0;
            }

            return childrenSize;
        }

        template <typename T>
        void Index<T>::SelectHeadDynamically(const std::shared_ptr<COMMON::BKTree> p_tree, int p_vectorCount, std::vector<int>& p_selected) {
            p_selected.clear();
            p_selected.reserve(p_vectorCount);

            if (static_cast<int>(std::round(m_options.m_ratio * p_vectorCount)) >= p_vectorCount)
            {
                for (int i = 0; i < p_vectorCount; ++i)
                {
                    p_selected.push_back(i);
                }

                return;
            }
            Options opts = m_options;

            int selectThreshold = m_options.m_selectThreshold;
            int splitThreshold = m_options.m_splitThreshold;

            double minDiff = 100;
            for (int select = 2; select <= m_options.m_selectThreshold; ++select)
            {
                opts.m_selectThreshold = select;
                opts.m_splitThreshold = m_options.m_splitThreshold;

                int l = m_options.m_splitFactor;
                int r = m_options.m_splitThreshold;

                while (l < r - 1)
                {
                    opts.m_splitThreshold = (l + r) / 2;
                    p_selected.clear();

                    SelectHeadDynamicallyInternal(p_tree, 0, opts, p_selected);
                    std::sort(p_selected.begin(), p_selected.end());
                    p_selected.erase(std::unique(p_selected.begin(), p_selected.end()), p_selected.end());

                    double diff = static_cast<double>(p_selected.size()) / p_vectorCount - m_options.m_ratio;

                    LOG(Helper::LogLevel::LL_Info,
                        "Select Threshold: %d, Split Threshold: %d, diff: %.2lf%%.\n",
                        opts.m_selectThreshold,
                        opts.m_splitThreshold,
                        diff * 100.0);

                    if (minDiff > fabs(diff))
                    {
                        minDiff = fabs(diff);

                        selectThreshold = opts.m_selectThreshold;
                        splitThreshold = opts.m_splitThreshold;
                    }

                    if (diff > 0)
                    {
                        l = (l + r) / 2;
                    }
                    else
                    {
                        r = (l + r) / 2;
                    }
                }
            }

            opts.m_selectThreshold = selectThreshold;
            opts.m_splitThreshold = splitThreshold;

            LOG(Helper::LogLevel::LL_Info,
                "Final Select Threshold: %d, Split Threshold: %d.\n",
                opts.m_selectThreshold,
                opts.m_splitThreshold);

            p_selected.clear();
            SelectHeadDynamicallyInternal(p_tree, 0, opts, p_selected);
            std::sort(p_selected.begin(), p_selected.end());
            p_selected.erase(std::unique(p_selected.begin(), p_selected.end()), p_selected.end());
        }

        template <typename T>
        bool Index<T>::SelectHead(std::shared_ptr<Helper::VectorSetReader>& p_reader) {
            std::shared_ptr<VectorSet> vectorset = p_reader->GetVectorSet();
            if (m_options.m_distCalcMethod == DistCalcMethod::Cosine && !p_reader->IsNormalized())
                vectorset->Normalize(m_options.m_iSelectHeadNumberOfThreads);
            COMMON::Dataset<T> data(vectorset->Count(), vectorset->Dimension(), vectorset->Count(), vectorset->Count() + 1, (T*)vectorset->GetData());

            auto t1 = std::chrono::high_resolution_clock::now();
            SelectHeadAdjustOptions(data.R());
            std::vector<int> selected;
            if (data.R() == 1) {
                selected.push_back(0);
            }
            else if (Helper::StrUtils::StrEqualIgnoreCase(m_options.m_selectType.c_str(), "Random")) {
                LOG(Helper::LogLevel::LL_Info, "Start generating Random head.\n");
                selected.resize(data.R());
                for (int i = 0; i < data.R(); i++) selected[i] = i;
                std::random_shuffle(selected.begin(), selected.end());
                int headCnt = static_cast<int>(std::round(m_options.m_ratio * data.R()));
                selected.resize(headCnt);
            }
            else if (Helper::StrUtils::StrEqualIgnoreCase(m_options.m_selectType.c_str(), "BKT")) {
                LOG(Helper::LogLevel::LL_Info, "Start generating BKT.\n");
                std::shared_ptr<COMMON::BKTree> bkt = std::make_shared<COMMON::BKTree>();
                bkt->m_iBKTKmeansK = m_options.m_iBKTKmeansK;
                bkt->m_iBKTLeafSize = m_options.m_iBKTLeafSize;
                bkt->m_iSamples = m_options.m_iSamples;
                bkt->m_iTreeNumber = m_options.m_iTreeNumber;
                bkt->m_fBalanceFactor = m_options.m_fBalanceFactor;
                LOG(Helper::LogLevel::LL_Info, "Start invoking BuildTrees.\n");
                LOG(Helper::LogLevel::LL_Info, "BKTKmeansK: %d, BKTLeafSize: %d, Samples: %d, BKTLambdaFactor:%f TreeNumber: %d, ThreadNum: %d.\n",
                    bkt->m_iBKTKmeansK, bkt->m_iBKTLeafSize, bkt->m_iSamples, bkt->m_fBalanceFactor, bkt->m_iTreeNumber, m_options.m_iSelectHeadNumberOfThreads);

                bkt->BuildTrees<T>(data, m_options.m_distCalcMethod, m_options.m_iSelectHeadNumberOfThreads, nullptr, nullptr, true);
                auto t2 = std::chrono::high_resolution_clock::now();
                double elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
                LOG(Helper::LogLevel::LL_Info, "End invoking BuildTrees.\n");
                LOG(Helper::LogLevel::LL_Info, "Invoking BuildTrees used time: %.2lf minutes (about %.2lf hours).\n", elapsedSeconds / 60.0, elapsedSeconds / 3600.0);

                if (m_options.m_saveBKT) {
                    std::stringstream bktFileNameBuilder;
                    bktFileNameBuilder << m_options.m_vectorPath << ".bkt." << m_options.m_iBKTKmeansK << "_"
                                       << m_options.m_iBKTLeafSize << "_" << m_options.m_iTreeNumber << "_" << m_options.m_iSamples << "_"
                                       << static_cast<int>(m_options.m_distCalcMethod) << ".bin";
                    bkt->SaveTrees(bktFileNameBuilder.str());
                }
                LOG(Helper::LogLevel::LL_Info, "Finish generating BKT.\n");

                LOG(Helper::LogLevel::LL_Info, "Start selecting nodes...Select Head Dynamically...\n");
                SelectHeadDynamically(bkt, data.R(), selected);

                if (selected.empty()) {
                    LOG(Helper::LogLevel::LL_Error, "Can't select any vector as head with current settings\n");
                    return false;
                }
            }

            LOG(Helper::LogLevel::LL_Info,
                "Seleted Nodes: %u, about %.2lf%% of total.\n",
                static_cast<unsigned int>(selected.size()),
                selected.size() * 100.0 / data.R());

            if (!m_options.m_noOutput)
            {
                std::sort(selected.begin(), selected.end());

                std::shared_ptr<Helper::DiskPriorityIO> output = SPTAG::f_createIO(), outputIDs = SPTAG::f_createIO();
                if (output == nullptr || outputIDs == nullptr ||
                    !output->Initialize((m_options.m_indexDirectory + FolderSep + m_options.m_headVectorFile).c_str(), std::ios::binary | std::ios::out) ||
                    !outputIDs->Initialize((m_options.m_indexDirectory + FolderSep + m_options.m_headIDFile).c_str(), std::ios::binary | std::ios::out)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to create output file:%s %s\n",
                        (m_options.m_indexDirectory + FolderSep + m_options.m_headVectorFile).c_str(),
                        (m_options.m_indexDirectory + FolderSep + m_options.m_headIDFile).c_str());
                    return false;
                }

                SizeType val = static_cast<SizeType>(selected.size());
                if (output->WriteBinary(sizeof(val), reinterpret_cast<char*>(&val)) != sizeof(val)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to write output file!\n");
                    return false;
                }
                DimensionType dt = data.C();
                if (output->WriteBinary(sizeof(dt), reinterpret_cast<char*>(&dt)) != sizeof(dt)) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to write output file!\n");
                    return false;
                }

                for (int i = 0; i < selected.size(); i++)
                {
                    uint64_t vid = static_cast<uint64_t>(selected[i]);
                    if (outputIDs->WriteBinary(sizeof(vid), reinterpret_cast<char*>(&vid)) != sizeof(vid)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to write output file!\n");
                        return false;
                    }

                    if (output->WriteBinary(sizeof(T) * data.C(), (char*)(data[vid])) != sizeof(T) * data.C()) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to write output file!\n");
                        return false;
                    }
                }
            }
            auto t3 = std::chrono::high_resolution_clock::now();
            double elapsedSeconds = std::chrono::duration_cast<std::chrono::seconds>(t3 - t1).count();
            LOG(Helper::LogLevel::LL_Info, "Total used time: %.2lf minutes (about %.2lf hours).\n", elapsedSeconds / 60.0, elapsedSeconds / 3600.0);
            return true;
        }

        template <typename T>
        ErrorCode Index<T>::BuildIndexInternal(std::shared_ptr<Helper::VectorSetReader>& p_reader) {
            if (!m_options.m_indexDirectory.empty()) {
                if (!direxists(m_options.m_indexDirectory.c_str()))
                {
                    mkdir(m_options.m_indexDirectory.c_str());
                }
            }

            auto t1 = std::chrono::high_resolution_clock::now();
            if (m_options.m_selectHead) {
                omp_set_num_threads(m_options.m_iSelectHeadNumberOfThreads);
                if (!SelectHead(p_reader)) {
                    LOG(Helper::LogLevel::LL_Error, "SelectHead Failed!\n");
                    return ErrorCode::Fail;
                }
            }
            auto t2 = std::chrono::high_resolution_clock::now();
            double selectHeadTime = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
            LOG(Helper::LogLevel::LL_Info, "select head time: %.2lfs\n", selectHeadTime);

            if (m_options.m_buildHead) {
                auto valueType = SPTAG::COMMON::DistanceUtils::Quantizer ? SPTAG::VectorValueType::UInt8 : m_options.m_valueType;
                m_index = SPTAG::VectorIndex::CreateInstance(m_options.m_indexAlgoType, valueType);
                m_index->SetParameter("DistCalcMethod", SPTAG::Helper::Convert::ConvertToString(m_options.m_distCalcMethod));
                for (const auto& iter : m_headParameters)
                {
                    m_index->SetParameter(iter.first.c_str(), iter.second.c_str());
                }

                std::shared_ptr<Helper::ReaderOptions> vectorOptions(new Helper::ReaderOptions(valueType, m_options.m_dim, VectorFileType::DEFAULT));
                auto vectorReader = Helper::VectorSetReader::CreateInstance(vectorOptions);
                if (ErrorCode::Success != vectorReader->LoadFile(m_options.m_indexDirectory + FolderSep + m_options.m_headVectorFile))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read head vector file.\n");
                    return ErrorCode::Fail;
                }
                if (m_index->BuildIndex(vectorReader->GetVectorSet(), nullptr, false, true) != ErrorCode::Success ||
                    m_index->SaveIndex(m_options.m_indexDirectory + FolderSep + m_options.m_headIndexFolder) != ErrorCode::Success) {
                    LOG(Helper::LogLevel::LL_Error, "Failed to build head index.\n");
                    return ErrorCode::Fail;
                }
            }
            auto t3 = std::chrono::high_resolution_clock::now();
            double buildHeadTime = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
            LOG(Helper::LogLevel::LL_Info, "select head time: %.2lfs build head time: %.2lfs\n", selectHeadTime, buildHeadTime);

            if (m_options.m_enableSSD) {
                omp_set_num_threads(m_options.m_iSSDNumberOfThreads);

                if (m_index == nullptr && LoadIndex(m_options.m_indexDirectory + FolderSep + m_options.m_headIndexFolder, m_index) != ErrorCode::Success) {
                    LOG(Helper::LogLevel::LL_Error, "Cannot load head index from %s!\n", (m_options.m_indexDirectory + FolderSep + m_options.m_headIndexFolder).c_str());
                    return ErrorCode::Fail;
                }
                if (!CheckHeadIndexType()) return ErrorCode::Fail;

                m_index->SetParameter("NumberOfThreads", std::to_string(m_options.m_iSSDNumberOfThreads));
                m_index->SetParameter("MaxCheck", std::to_string(m_options.m_maxCheck));
                m_index->SetParameter("HashTableExponent", std::to_string(m_options.m_hashExp));
                m_index->UpdateIndex();

                if (m_options.m_useKV)
                {
                    m_extraSearcher.reset(new ExtraRocksDBController<T>(m_options.m_KVPath.c_str()));
                } else {
                    m_extraSearcher.reset(new ExtraFullGraphSearcher<T>());
                }
                if (m_options.m_buildSsdIndex) {
                    if (!m_extraSearcher->BuildIndex(p_reader, m_index, m_options)) {
                        LOG(Helper::LogLevel::LL_Error, "BuildSSDIndex Failed!\n");
                        return ErrorCode::Fail;
                    }
                }
                if (!m_extraSearcher->LoadIndex(m_options)) {
                    LOG(Helper::LogLevel::LL_Error, "Cannot Load SSDIndex!\n");
                    return ErrorCode::Fail;
                }

                if (!m_options.m_useKV) {
                    m_vectorTranslateMap.reset(new std::uint64_t[m_index->GetNumSamples()], std::default_delete<std::uint64_t[]>());
                    std::shared_ptr<Helper::DiskPriorityIO> ptr = SPTAG::f_createIO();
                    if (ptr == nullptr || !ptr->Initialize((m_options.m_indexDirectory + FolderSep + m_options.m_headIDFile).c_str(), std::ios::binary | std::ios::in)) {
                        LOG(Helper::LogLevel::LL_Error, "Failed to open headIDFile file:%s\n", (m_options.m_indexDirectory + FolderSep + m_options.m_headIDFile).c_str());
                        return ErrorCode::Fail;
                    }
                    IOBINARY(ptr, ReadBinary, sizeof(std::uint64_t) * m_index->GetNumSamples(), (char*)(m_vectorTranslateMap.get()));
                } else {
                    //data structrue initialization
                    m_deletedID.Load(m_options.m_fullDeletedIDFile, m_iDataBlockSize, m_iDataCapacity);
                    m_rwLocks = std::make_unique<std::shared_timed_mutex[]>(500000000);
                    m_postingSizes = std::make_unique<std::atomic_uint32_t[]>(500000000);
                    std::ifstream input(m_options.m_ssdInfoFile, std::ios::binary);
                    if (!input.is_open())
                    {
                        fprintf(stderr, "Failed to open file: %s\n", m_options.m_ssdInfoFile.c_str());
                        exit(1);
                    }

					int vectorNum;
                    input.read(reinterpret_cast<char*>(&vectorNum), sizeof(vectorNum));
					m_vectorNum.store(vectorNum);

					LOG(Helper::LogLevel::LL_Info, "Current vector num: %d.\n", m_vectorNum.load());

					uint32_t postingNum;
					input.read(reinterpret_cast<char*>(&postingNum), sizeof(postingNum));

					LOG(Helper::LogLevel::LL_Info, "Current posting num: %d.\n", postingNum);

					for (int idx = 0; idx < postingNum; idx++) {
						uint32_t tmp;
						input.read(reinterpret_cast<char*>(&tmp), sizeof(uint32_t));
						m_postingSizes[idx].store(tmp);
					}

					input.close();    
                }       
            }
            auto t4 = std::chrono::high_resolution_clock::now();
            double buildSSDTime = std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count();
            LOG(Helper::LogLevel::LL_Info, "select head time: %.2lfs build head time: %.2lfs build ssd time: %.2lfs\n", selectHeadTime, buildHeadTime, buildSSDTime);

            if (m_options.m_deleteHeadVectors) {
                if (fileexists((m_options.m_indexDirectory + FolderSep + m_options.m_headVectorFile).c_str()) &&
                    remove((m_options.m_indexDirectory + FolderSep + m_options.m_headVectorFile).c_str()) != 0) {
                    LOG(Helper::LogLevel::LL_Warning, "Head vector file can't be removed.\n");
                }
            }

            m_workSpacePool.reset(new COMMON::WorkSpacePool<ExtraWorkSpace>());
            m_workSpacePool->Init(m_options.m_iSSDNumberOfThreads, m_options.m_maxCheck, m_options.m_hashExp, m_options.m_searchInternalResultNum, min(m_options.m_postingPageLimit, m_options.m_searchPostingPageLimit + 1) << PageSizeEx);
            m_bReady = true;
            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::BuildIndex(bool p_normalized)
        {
            SPTAG::VectorValueType valueType = SPTAG::COMMON::DistanceUtils::Quantizer ? SPTAG::VectorValueType::UInt8 : m_options.m_valueType;
            std::shared_ptr<Helper::ReaderOptions> vectorOptions(new Helper::ReaderOptions(valueType, m_options.m_dim, m_options.m_vectorType, m_options.m_vectorDelimiter, p_normalized));
            auto vectorReader = Helper::VectorSetReader::CreateInstance(vectorOptions);
            if (m_options.m_vectorPath.empty())
            {
                LOG(Helper::LogLevel::LL_Info, "Vector file is empty. Skipping loading.\n");
            }
            else {
                if (ErrorCode::Success != vectorReader->LoadFile(m_options.m_vectorPath))
                {
                    LOG(Helper::LogLevel::LL_Error, "Failed to read vector file.\n");
                    return ErrorCode::Fail;
                }
                m_options.m_vectorSize = vectorReader->GetVectorSet()->Count();
            }

            return BuildIndexInternal(vectorReader);
        }

        template <typename T>
        ErrorCode Index<T>::BuildIndex(const void* p_data, SizeType p_vectorNum, DimensionType p_dimension, bool p_normalized)
        {
            if (p_data == nullptr || p_vectorNum == 0 || p_dimension == 0) return ErrorCode::EmptyData;

            if (m_options.m_distCalcMethod == DistCalcMethod::Cosine && !p_normalized) {
                COMMON::Utils::BatchNormalize((T*)p_data, p_vectorNum, p_dimension, COMMON::Utils::GetBase<T>(), m_options.m_iSSDNumberOfThreads);
            }
            std::shared_ptr<VectorSet> vectorSet(new BasicVectorSet(ByteArray((std::uint8_t*)p_data, p_vectorNum * p_dimension * sizeof(T), false),
                                                                    GetEnumValueType<T>(), p_dimension, p_vectorNum));
            SPTAG::VectorValueType valueType = SPTAG::COMMON::DistanceUtils::Quantizer ? SPTAG::VectorValueType::UInt8 : m_options.m_valueType;
            std::shared_ptr<Helper::VectorSetReader> vectorReader(new Helper::MemoryVectorReader(std::make_shared<Helper::ReaderOptions>(valueType, p_dimension, VectorFileType::DEFAULT, m_options.m_vectorDelimiter, m_options.m_iSSDNumberOfThreads, true),
                                                                                                 vectorSet));

            m_options.m_vectorSize = p_vectorNum;
            return BuildIndexInternal(vectorReader);
        }

        template <typename T>
        ErrorCode Index<T>::UpdateIndex()
        {
            omp_set_num_threads(m_options.m_iSSDNumberOfThreads);
            m_index->UpdateIndex();
            m_workSpacePool.reset(new COMMON::WorkSpacePool<ExtraWorkSpace>());
            m_workSpacePool->Init(m_options.m_iSSDNumberOfThreads, m_options.m_maxCheck, m_options.m_hashExp, m_options.m_searchInternalResultNum, min(m_options.m_postingPageLimit, m_options.m_searchPostingPageLimit + 1) << PageSizeEx);
            return ErrorCode::Success;
        }

        template <typename T>
        ErrorCode Index<T>::SetParameter(const char* p_param, const char* p_value, const char* p_section)
        {
            if (SPTAG::Helper::StrUtils::StrEqualIgnoreCase(p_section, "BuildHead") && !SPTAG::Helper::StrUtils::StrEqualIgnoreCase(p_param, "isExecute")) {
                if (m_index != nullptr) return m_index->SetParameter(p_param, p_value);
                else m_headParameters[p_param] = p_value;
            }
            else {
                m_options.SetParameter(p_section, p_param, p_value);
            }
            if (SPTAG::Helper::StrUtils::StrEqualIgnoreCase(p_param, "DistCalcMethod")) {
                m_fComputeDistance = COMMON::DistanceCalcSelector<T>(m_options.m_distCalcMethod);
                m_iBaseSquare = (m_options.m_distCalcMethod == DistCalcMethod::Cosine) ? COMMON::Utils::GetBase<T>() * COMMON::Utils::GetBase<T>() : 1;
            }
            return ErrorCode::Success;
        }

        template <typename T>
        std::string Index<T>::GetParameter(const char* p_param, const char* p_section) const
        {
            if (SPTAG::Helper::StrUtils::StrEqualIgnoreCase(p_section, "BuildHead") && !SPTAG::Helper::StrUtils::StrEqualIgnoreCase(p_param, "isExecute")) {
                if (m_index != nullptr) return m_index->GetParameter(p_param);
                else {
                    auto iter = m_headParameters.find(p_param);
                    if (iter != m_headParameters.end()) return iter->second;
                    return "Undefined!";
                }
            }
            else {
                return m_options.GetParameter(p_section, p_param);
            }
        }

        // Add insert entry to persistent buffer
        template <typename T>
        ErrorCode Index<T>::AddIndex(const void *p_data, SizeType p_vectorNum, DimensionType p_dimension,
                                     std::shared_ptr<MetadataSet> p_metadataSet, bool p_withMetaIndex,
                                     bool p_normalized)
        {
            if (m_options.m_indexAlgoType != IndexAlgoType::BKT || m_extraSearcher == nullptr) {
                LOG(Helper::LogLevel::LL_Error, "Only Support BKT Update");
                return ErrorCode::Fail;
            }

            std::vector<QueryResult> p_queryResults(p_vectorNum, QueryResult(NULL, m_options.m_internalResultNum, false));

            for (int k = 0; k < p_vectorNum; k++)
            {
                p_queryResults[k].SetTarget(reinterpret_cast<const T*>(reinterpret_cast<const char*>(p_data) + k * p_dimension));
                p_queryResults[k].Reset();
                auto VID = m_vectorNum++;
                {
                    std::lock_guard<std::mutex> lock(m_dataAddLock);
                    m_deletedID.AddBatch(1);
                    m_reassignedID.AddBatch(1);
                }

                m_index->SearchIndex(p_queryResults[k]);

                int replicaCount = 0;
                BasicResult* queryResults = p_queryResults[k].GetResults();
                std::vector<EdgeInsert> selections(static_cast<size_t>(m_options.m_replicaCount));
                for (int i = 0; i < p_queryResults[k].GetResultNum() && replicaCount < m_options.m_replicaCount; ++i)
                {
                    if (queryResults[i].VID == -1) {
                        break;
                    }
                    // RNG Check.
                    bool rngAccpeted = true;
                    for (int j = 0; j < replicaCount; ++j)
                    {
                        float nnDist = m_index->ComputeDistance(m_index->GetSample(queryResults[i].VID),
                                                                m_index->GetSample(selections[j].headID));
                        if (nnDist <= queryResults[i].Dist)
                        {
                            rngAccpeted = false;
                            break;
                        }
                    }
                    if (!rngAccpeted)
                        continue;
                    selections[replicaCount].headID = queryResults[i].VID;
                    selections[replicaCount].fullID = VID;
                    selections[replicaCount].distance = queryResults[i].Dist;
                    selections[replicaCount].order = (char)replicaCount;
                    ++replicaCount;
                }

                char insertCode = 0;
                SizeType assignID = 0;
                for (int i = 0; i < replicaCount; i++)
                {
                    std::string assignment;
                    assignment += Helper::Convert::Serialize<char>(&insertCode, 1);
                    assignment += Helper::Convert::Serialize<int>(&selections[i].headID, 1);
                    assignment += Helper::Convert::Serialize<int>(&VID, 1);
                    assignment += Helper::Convert::Serialize<T>(p_queryResults[k].GetTarget(), m_options.m_dim);
                    assignID = m_persistentBuffer->PutAssignment(assignment);
                }
            }

        }

        // Add delete entry to persistent buffer
        template <typename T>
        ErrorCode Index<T>::DeleteIndex(const SizeType &p_id)
        {
            // I think we should just delete instantly
            if (m_options.m_addDeleteTaskToPM) {
                char deleteCode = 1;
                int VID = p_id;
                std::string assignment;
                assignment += Helper::Convert::Serialize<char>(&deleteCode, 1);
                assignment += Helper::Convert::Serialize<int>(&VID, 1);
                m_persistentBuffer->PutAssignment(assignment);
            } else {
                std::lock_guard<std::mutex> lock(m_dataAddLock);
                m_deletedID.Insert(p_id);
            }
            return ErrorCode::Success;
        }

        template <typename T>
        void SPTAG::SPANN::Index<T>::Dispatcher::dispatch()
        {
            while (running) {
                bool noAssignment = true;
                int currentAssignmentID = m_persistentBuffer->GetCurrentAssignmentID();
                int scanNum = std::min<int>(sentAssignment + batch, currentAssignmentID);
                if (scanNum != sentAssignment) {
                    noAssignment = false;
                }

                std::map<SizeType, std::shared_ptr<std::string>> newPart;
                newPart.clear();
                int i;
                for (i = sentAssignment; i < scanNum; i++) {
                    std::string assignment;
                    m_persistentBuffer->GetAssignment(i, &assignment);
                    if(assignment.size() == 0)
                        break;
                    uint8_t* postingP = reinterpret_cast<uint8_t*>(&assignment.front());
                    char code = *(reinterpret_cast<char*>(postingP));
                    if (assignment.size() == 0)
                    {
                        uint8_t* headPointer = postingP + sizeof(char);
                        LOG(Helper::LogLevel::LL_Info, "Error\n");
                        LOG(Helper::LogLevel::LL_Info, "code: %d, headID: %d, assignment size: %d\n", code, *(reinterpret_cast<int*>(headPointer)), assignment.size());
                        LOG(Helper::LogLevel::LL_Info, "ScanNum: %d, SentNum: %d, CurrentAssignNum: %d, ProcessingAssignment: %d\n", scanNum, sentAssignment.load(), currentAssignmentID, i);
                        exit(1);
                    }
                    if (code == 0) {
                        // insert
                        uint8_t* headPointer = postingP + sizeof(char);
                        int32_t headID = *(reinterpret_cast<int*>(headPointer));
                        int32_t vid = *(reinterpret_cast<int*>(headPointer + sizeof(int)));
                        if (m_index->CheckIdDeleted(vid)) {
                            continue;
                        }
                        if (newPart.find(headID) == newPart.end()) {
                            newPart[headID].reset(new std::string(Helper::Convert::Serialize<uint8_t>(headPointer + sizeof(int), m_index->GetValueSize() + sizeof(int))));
                        } else {
                            *newPart[headID] += Helper::Convert::Serialize<uint8_t>(headPointer + sizeof(int), m_index->GetValueSize() + sizeof(int));
                        }
                    } else {
                        // delete
                        uint8_t* vectorPointer = postingP + sizeof(char);
                        int VID = *(reinterpret_cast<int*>(vectorPointer));
                        //LOG(Helper::LogLevel::LL_Info, "Scanner: delete: %d\n", VID);
                        m_index->DeleteIndex(VID);
                    }
                }

                for (auto iter = newPart.begin(); iter != newPart.end(); iter++) {
                    int appendNum = (*iter->second).size() / (m_index->GetValueSize() + sizeof(int));
                    m_index->AppendAsync(iter->first, appendNum, iter->second);
                }

                sentAssignment = i;
                if (noAssignment) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                } else {
                    //LOG(Helper::LogLevel::LL_Info, "Process Append Assignments: %d, Delete Assignments: %d\n", newPart.size(), deletedVector.size());
                }
            }
        }

        template <typename ValueType>
        ErrorCode SPTAG::SPANN::Index<ValueType>::Split(const SizeType headID, int appendNum, std::string& appendPosting)
        {
            // TimeUtils::StopW sw;
            std::unique_lock<std::shared_timed_mutex> lock(m_rwLocks[headID]);
            if (m_postingSizes[headID].load() + appendNum < m_options.m_postingVectorLimit) {
                return ErrorCode::FailSplit;
            }
            m_splitTaskNum++;
            std::string postingList;
            m_extraSearcher->SearchIndex(headID, postingList);
            postingList += appendPosting;

            // reinterpret postingList to vectors and IDs
            auto* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
            int postVectorNum = postingList.size() / (m_options.m_vectorSize + sizeof(int));
            COMMON::Dataset<ValueType> smallSample;  // smallSample[i] -> VID
            std::shared_ptr<uint8_t> vectorBuffer(new uint8_t[m_options.m_vectorSize * postVectorNum], std::default_delete<uint8_t[]>());
            std::vector<int> localIndicesInsert(postVectorNum);  // smallSample[i] = j <-> localindices[j] = i
            std::vector<int> localIndices(postVectorNum);
            auto vectorBuf = vectorBuffer.get();
            int realVectorNum = postVectorNum;
            int index = 0;
            //LOG(Helper::LogLevel::LL_Info, "Scanning\n");
            for (int j = 0; j < postVectorNum; j++)
            {
                uint8_t* vectorId = postingP + j * (m_options.m_vectorSize + sizeof(int));
                //LOG(Helper::LogLevel::LL_Info, "vector index/total:id: %d/%d:%d\n", j, m_postingSizes[selections[i].headID], *(reinterpret_cast<int*>(vectorId)));
                if (CheckIdDeleted(*(reinterpret_cast<int*>(vectorId)))) {
                    realVectorNum--;
                } else {
                    localIndicesInsert[index] = *(reinterpret_cast<int*>(vectorId));
                    localIndices[index] = index;
                    index++;
                    memcpy(vectorBuf, vectorId + sizeof(int), m_options.m_vectorSize);
                    vectorBuf += m_options.m_vectorSize;
                }
            }
            // double gcEndTime = sw.getElapsedMs();
            // m_splitGcCost += gcEndTime;
            if (realVectorNum < m_options.m_postingVectorLimit)
            {
                postingList.clear();
                for (int j = 0; j < realVectorNum; j++)
                {
                    postingList += Helper::Convert::Serialize<int>(&localIndicesInsert[j], 1);
                    postingList += Helper::Convert::Serialize<ValueType>(vectorBuffer.get() + j * m_options.m_vectorSize, m_options.m_dim);
                }
                m_postingSizes[headID].store(realVectorNum);
                m_extraSearcher->AddIndex(headID, postingList);
                // m_splitWriteBackCost += sw.getElapsedMs() - gcEndTime;
                return ErrorCode::Success;
            }
            //LOG(Helper::LogLevel::LL_Info, "Resize\n");
            localIndicesInsert.resize(realVectorNum);
            localIndices.resize(realVectorNum);
            smallSample.Initialize(realVectorNum, m_options.m_dim, m_iDataBlockSize, m_iDataCapacity, reinterpret_cast<ValueType*>(vectorBuffer.get()), false);

            //LOG(Helper::LogLevel::LL_Info, "Headid: %d Sample Vector Num: %d, Real Vector Num: %d\n", selections[i].headID, smallSample.R(), realVectorNum);

            // k = 2, maybe we can change the split number, now it is fixed
            SPTAG::COMMON::KmeansArgs<ValueType> args(2, smallSample.C(), (SizeType)localIndicesInsert.size(), 1, m_index->GetDistCalcMethod());
            std::shuffle(localIndices.begin(), localIndices.end(), std::mt19937(std::random_device()()));
            int numClusters = SPTAG::COMMON::KmeansClustering(smallSample, localIndices, 0, (SizeType)localIndices.size(), args);
            if (numClusters <= 1)
            {
                postingList.clear();
                float r = 0.f;
                for (int j = 0; j < realVectorNum; j++)
                {
                    postingList += Helper::Convert::Serialize<int>(&localIndicesInsert[j], 1);
                    postingList += Helper::Convert::Serialize<ValueType>(vectorBuffer.get() + j * m_options.m_vectorSize, m_options.m_dim);
                    auto dist = m_index->ComputeDistance(vectorBuffer.get() + j * m_options.m_vectorSize, m_index->GetSample(headID));
                    r = std::max<float>(r, dist);
                }
                m_postingSizes[headID].store(realVectorNum);
                m_extraSearcher->AddIndex(headID, postingList);
                return ErrorCode::Success;
            }
            // double clusterEndTime = sw.getElapsedMs();
            // m_splitClusteringCost += clusterEndTime - gcEndTime;

            long long newHeadVID = -1;
            int first = 0;
            std::vector<SizeType> newHeadsID(2);
            std::vector<std::string> newPostingLists;
            for (int k = 0; k < 2; k++) {
                int begin, end = 0;
                std::string postingList;
                if (args.counts[k] == 0)	continue;

                // LOG(Helper::LogLevel::LL_Info, "Insert new head vector\n");
                // Notice: newHeadVID maybe an existing head vector

                m_index->AddIndexId(smallSample[args.clusterIdx[k]], 1, m_options.m_dim, begin, end);
                newHeadVID = begin;
                newHeadsID.push_back(begin);
                //LOG(Helper::LogLevel::LL_Info, "Head id: %d split into : %d\n", headID, newHeadVID);

                for (int j = 0; j < args.counts[k]; j++)
                {
                    postingList += Helper::Convert::Serialize<SizeType>(&localIndicesInsert[localIndices[first + j]], 1);
                    postingList += Helper::Convert::Serialize<ValueType>(smallSample[localIndices[first + j]], m_options.m_dim);
                }
                m_extraSearcher->AddIndex(newHeadVID, postingList);
                newPostingLists.push_back(postingList);
                first += args.counts[k];

                // TODO: move postingSizes to extraSearcher
                m_postingSizes[newHeadVID] = args.counts[k];
                m_index->AddIndexIdx(begin, end);
            }

            m_index->DeleteIndex(headID);
            m_extraSearcher->DeleteIndex(headID);

            m_postingSizes[headID] = 0;
            lock.unlock();
            // int split_order = ++m_split_num;
            // m_splitUpdateIndexCost += sw.getElapsedMs() - clusterEndTime;
            //QuantifySplit(headID, newPostingLists, newHeadsID, split_order);

            if (!m_options.m_disableReassign) ReAssign(headID, newPostingLists, newHeadsID);

            //LOG(Helper::LogLevel::LL_Info, "After ReAssign\n");

            //QuantifySplit(headID, newPostingLists, newHeadsID, split_order);
            return ErrorCode::Success;
        }

        template <typename ValueType>
        ErrorCode SPTAG::SPANN::Index<ValueType>::ReAssign(SizeType headID, std::vector<std::string>& postingLists, std::vector<SizeType>& newHeadsID) {
//            TimeUtils::StopW sw;
            auto headVector = reinterpret_cast<const ValueType*>(m_index->GetSample(headID));
            if (m_options.m_reassignK > 0) {
                COMMON::QueryResultSet<ValueType> nearbyHeads(NULL, m_options.m_reassignK);
                nearbyHeads.SetTarget(headVector);
                nearbyHeads.Reset();
                m_index->SearchIndex(nearbyHeads);
                BasicResult* queryResults = nearbyHeads.GetResults();
                postingLists.resize(nearbyHeads.GetResultNum() + postingLists.size());
                for (int i = 0; i < nearbyHeads.GetResultNum(); i++) {
                    auto vid = queryResults[i].VID;
                    if (vid == -1) {
                        break;
                    }
                    if (find(newHeadsID.begin(), newHeadsID.end(), vid) == newHeadsID.end()) {
                        m_extraSearcher->SearchIndex(vid, postingLists[i + 2]);
                    }
                }
//                m_reassignSearchHeadCost += sw.getElapsedMs();
            }

            std::map<SizeType, ValueType*> reAssignVectorsTop0;
            std::map<SizeType, ValueType*> reAssignVectorsTopK;
            for (int i = 0; i < postingLists.size(); i++) {
                auto& postingList = postingLists[i];
                int postVectorNum = postingList.size() / (m_options.m_vectorSize + sizeof(int));
                auto* postingP = reinterpret_cast<uint8_t*>(&postingList.front());
                for (int j = 0; j < postVectorNum; j++) {
                    uint8_t* vectorId = postingP + j * (m_options.m_vectorSize + sizeof(int));
                    SizeType vid = *(reinterpret_cast<SizeType*>(vectorId));
                    if (i <= 1) {
                        if (reAssignVectorsTop0.find(vid) == reAssignVectorsTop0.end() && !CheckIdDeleted(vid))
                            reAssignVectorsTop0[vid] = reinterpret_cast<ValueType*>(vectorId + sizeof(int));
                        //PrintFirstFiveDimInt8(vectorId + sizeof(int), vid);
                    } else {
                        if (reAssignVectorsTop0.find(vid) == reAssignVectorsTop0.end())
                        {
                            if (reAssignVectorsTopK.find(vid) == reAssignVectorsTopK.end() && !CheckIdDeleted(vid))
                                reAssignVectorsTopK[vid] = reinterpret_cast<ValueType*>(vectorId + sizeof(int));
                        }
                    }
                }
            }

            ReAssignVectors(reAssignVectorsTop0, newHeadsID);
            ReAssignVectors(reAssignVectorsTopK, newHeadsID, true);
//            m_reassignTotalCost += sw.getElapsedMs();
            return ErrorCode::Success;
        }

        template <typename ValueType>
        void SPTAG::SPANN::Index<ValueType>::ReAssignVectors(std::map<SizeType, ValueType*>& reAssignVectors,
                             std::vector<SizeType>& newHeadsID, bool check)
        {
            auto numQueries = reAssignVectors.size();
            if (!check) {
                {
                    std::lock_guard<std::mutex> lock(m_dataAddLock);
                    m_deletedID.AddBatch(numQueries);

                    m_reassignedID.AddBatch(numQueries);
                }
                int count = 0;
                auto newFirstVID = m_vectorNum.fetch_add(numQueries);
                for (auto it = reAssignVectors.begin(); it != reAssignVectors.end(); ++it) {
                    //m_currerntReassignTaskNum++;
                    //PrintFirstFiveDimInt8(reinterpret_cast<uint8_t*>(it->second), it->first);
                    std::shared_ptr<std::string> vectorContain(new std::string(Helper::Convert::Serialize<uint8_t>(it->second, m_options.m_vectorSize)));
                    //PrintFirstFiveDimInt8(reinterpret_cast<uint8_t*>(&vectorContain->front()), it->first);
                    ReassignAsync(vectorContain, newFirstVID + count, newHeadsID, false, it->first);
                    count++;
                }
            } else {
                for (auto it = reAssignVectors.begin(); it != reAssignVectors.end(); ++it) {
                    std::shared_ptr<std::string> vectorContain(new std::string(Helper::Convert::Serialize<uint8_t>(it->second, m_options.m_vectorSize)));
                    ReassignAsync(vectorContain, 0, newHeadsID, true, it->first);
                }
            }

//            m_reassigned += numQueries;
        }

        template <typename ValueType>
        void SPTAG::SPANN::Index<ValueType>::ReAssignUpdate
                (std::shared_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType> &newHeads, bool check,
                 SizeType oldVID)
        {
//            TimeUtils::StopW sw;

            COMMON::QueryResultSet<ValueType> p_queryResults(NULL, m_options.m_internalResultNum);
            //PrintFirstFiveDimInt8(reinterpret_cast<uint8_t*>(&vectorContain->front()), oldVID);
            p_queryResults.SetTarget(reinterpret_cast<ValueType*>(&vectorContain->front()));
            p_queryResults.Reset();
            m_index->SearchIndex(p_queryResults);

//            double indexSearchEndTime = sw.getElapsedMs();
//            m_reassignSearchVecotrCost += indexSearchEndTime;

            int replicaCount = 0;
            BasicResult* queryResults = p_queryResults.GetResults();
            std::vector<EdgeInsert> selections(static_cast<size_t>(m_options.m_replicaCount));
            if (check) {
                bool reassign = false;
                for (int i = 0; i < p_queryResults.GetResultNum(); i++) {
                    auto vid = queryResults[i].VID;
                    if (find(newHeads.begin(), newHeads.end(), vid) != newHeads.end()) {
                        reassign = true;
                        break;
                    } else if (vid == -1) {
                        break;
                    }
                }
                if (!reassign) {
                    return;
                }
                VID = m_vectorNum.fetch_add(1);
                {
                    std::lock_guard<std::mutex> lock(m_dataAddLock);
                    m_deletedID.AddBatch(1);
                    m_reassignedID.AddBatch(1);
                }
            }

            int i;
            for (i = 0; i < p_queryResults.GetResultNum() && replicaCount < m_options.m_replicaCount; ++i) {
                if (queryResults[i].VID == -1) {
                    break;
                }
                // RNG Check.
                bool rngAccpeted = true;
                for (int j = 0; j < replicaCount; ++j) {
                    float nnDist = m_index->ComputeDistance(
                            m_index->GetSample(queryResults[i].VID),
                            m_index->GetSample(selections[j].headID));
                    if (nnDist <= queryResults[i].Dist) {
                        rngAccpeted = false;
                        break;
                    }
                }
                if (!rngAccpeted)
                    continue;

                selections[replicaCount].headID = queryResults[i].VID;
                selections[replicaCount].fullID = VID;
                selections[replicaCount].distance = queryResults[i].Dist;
                selections[replicaCount].order = (char)replicaCount;
                ++replicaCount;
            }

//            double selectionEndTime = sw.getElapsedMs();
//            m_reassignSelectionCost += selectionEndTime - indexSearchEndTime;

            //LOG(Helper::LogLevel::LL_Info, "Reassign: oldVID:%d, replicaCount:%d, candidateNum:%d, dist0:%f\n", oldVID, replicaCount, i, selections[0].distance);

            for (i = 0; i < replicaCount; i++) {
                if (CheckIdDeleted(VID)) {
                    break;
                }
                std::shared_ptr<std::string> newPart(new std::string);
                *newPart += Helper::Convert::Serialize<int>(&VID, 1);
                *newPart += Helper::Convert::Serialize<ValueType>(p_queryResults.GetTarget(), m_options.m_dim);
                auto headID = selections[i].headID;
                //LOG(Helper::LogLevel::LL_Info, "Reassign: headID :%d, oldVID:%d, newVID:%d, posting length: %d, dist: %f\n", headID, oldVID, VID, m_postingSizes[headID].load(), selections[i].distance);
                Append(headID, 1, newPart, oldVID);
            }
//            m_reassignSsdCost += sw.getElapsedMs() - selectionEndTime;
            m_deletedID.Insert(oldVID);
        }

        template <typename ValueType>
        ErrorCode SPTAG::SPANN::Index<ValueType>::Append(SizeType headID, int appendNum, std::shared_ptr<std::string> appendPosting, SizeType oldVID)
        {
//            TimeUtils::StopW sw;
            m_appendTaskNum++;
            if (appendNum == 0) {
                LOG(Helper::LogLevel::LL_Info, "Error!, headID :%d, appendNum:%d\n", headID, appendNum);
            }
        checkDeleted:
            if (!m_index->ContainSample(headID)) {
                m_headMiss++;
                int newVID = m_vectorNum.fetch_add(appendNum);
                {
                    std::lock_guard<std::mutex> lock(m_dataAddLock);
                    m_deletedID.AddBatch(appendNum);
                    m_reassignedID.AddBatch(appendNum);
                }
                uint8_t* postingP = reinterpret_cast<uint8_t*>(&appendPosting->front());
                std::vector<SizeType> newHeads;
                for (int i = 0; i < appendNum; i++)
                {
//                    m_currerntReassignTaskNum++;
                    uint8_t* vid = postingP +  i * (m_options.m_vectorSize + sizeof(int));
                    std::shared_ptr<std::string> vectorContain(new std::string(Helper::Convert::Serialize<uint8_t>(vid + sizeof(int), m_options.m_vectorSize)));
                    ReassignAsync(vectorContain, newVID + i, newHeads, false, *(reinterpret_cast<int*>(vid)));
                }
                return ErrorCode::Success;
            }
            if (m_postingSizes[headID].load() + appendNum > m_options.m_postingVectorLimit) {
                // double splitStartTime = sw.getElapsedMs();
                if (Split(headID, appendNum, *appendPosting) == ErrorCode::FailSplit) {
                    goto checkDeleted;
                }
                // m_splitTotalCost += sw.getElapsedMs() - splitStartTime;
            } else {
                // double appendSsdStartTime = sw.getElapsedMs();
                {
                    std::shared_lock<std::shared_timed_mutex> lock(m_rwLocks[headID]);
                    if (!m_index->ContainSample(headID)) {
                        goto checkDeleted;
                    }
                    //LOG(Helper::LogLevel::LL_Info, "Merge: headID: %d, appendNum:%d\n", headID, appendNum);
                    m_extraSearcher->AppendPosting(headID, *appendPosting);
                }
                m_postingSizes[headID].fetch_add(appendNum, std::memory_order_relaxed);
                // m_appendSsdCost += sw.getElapsedMs() - appendSsdStartTime;
            }

            // m_appendTotalCost += sw.getElapsedMs();
            return ErrorCode::Success;
        }

        template <typename T>
        void SPTAG::SPANN::Index<T>::ProcessAsyncReassign(std::shared_ptr<std::string> vectorContain, SizeType VID, std::vector<SizeType>& newHeads, bool check,
                                                          SizeType oldVID, std::function<void()> p_callback)
        {
            //LOG(Helper::LogLevel::LL_Info, "ReassignID: %d, newID: %d\n", oldVID, VID);

            if (m_reassignedID.Contains(oldVID)) {
                m_deletedID.Insert(VID);
                return;
            }

            m_reassignedID.Insert(oldVID);

            ReAssignUpdate(vectorContain, VID, newHeads, check, oldVID);

            if (p_callback != nullptr) {
                p_callback();
            }
        }
    }
}

#define DefineVectorValueType(Name, Type) \
template class SPTAG::SPANN::Index<Type>; \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType


