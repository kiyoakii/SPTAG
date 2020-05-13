#ifdef DefineSSDParameter

// Both Building and Searching
DefineSSDParameter(m_buildSsdIndex, bool, true, "BuildSsdIndex")
DefineSSDParameter(m_vectorIDTranslate, string, string(""), "VectorIDTranslate")
DefineSSDParameter(m_headIndexFolder, string, string(""), "HeadIndexFolder")
DefineSSDParameter(m_queryFile, string, string(""), "QueryFile")
DefineSSDParameter(m_internalResultNum, int, 0, "InternalResultNum")
DefineSSDParameter(m_iNumberOfThreads, int, 16L, "NumberOfThreads")
DefineSSDParameter(m_headConfig, string, string(""), "HeadConfig") // must be in "Index" section
// query file meta
DefineSSDParameter(m_queryFileType, SPTAG::VectorFileType, SPTAG::VectorFileType::Undefined, "QueryFileType")
DefineSSDParameter(m_iQueryNumber, SPTAG::SizeType, -1, "QueryNumber")
DefineSSDParameter(m_iQueryDimension, SPTAG::DimensionType, -1, "QueryDimension")


// Building
DefineSSDParameter(m_ssdIndex, string, string(""), "SsdIndex")
DefineSSDParameter(m_replicaCount, int, 4, "ReplicaCount")
DefineSSDParameter(m_postingPageLimit, int, 2, "PostingPageLimit")
DefineSSDParameter(m_outputEmptyReplicaID, bool, false, "OutputEmptyReplicaID")


// Searching
DefineSSDParameter(m_searchResult, string, string(""), "SearchResult")
DefineSSDParameter(m_extraFullGraphFile, string, string(""), "ExtraFullGraphFile")
DefineSSDParameter(m_truthFile, string, string(""), "TruthFile")
DefineSSDParameter(m_warmupFile, string, string(""), "WarmupFile")
DefineSSDParameter(m_logFile, string, string(""), "LogFile")
DefineSSDParameter(m_extraMaxCheck, string, string("1024"), "ExtraMaxCheck")
DefineSSDParameter(m_qpsLimit, int, 0, "QpsLimit")
DefineSSDParameter(m_resultNum, int, 32, "ResultNum")
DefineSSDParameter(m_queryCountLimit, int, (numeric_limits<int>::max)(), "QueryCountLimit")
// warmup file meta
DefineSSDParameter(m_warmupFileType, SPTAG::VectorFileType, SPTAG::VectorFileType::Undefined, "WarmupFileType")
DefineSSDParameter(m_iWarmupNumber, SPTAG::SizeType, -1, "WarmupNumber")
DefineSSDParameter(m_iWarmupDimension, SPTAG::DimensionType, -1, "WarmupDimension")
// truth file meta
DefineSSDParameter(m_truthFileType, SPTAG::TruthFileType, SPTAG::TruthFileType::Undefined, "TruthFileType")
DefineSSDParameter(m_iTruthNumber, SPTAG::SizeType, -1, "TruthNumber")

#endif