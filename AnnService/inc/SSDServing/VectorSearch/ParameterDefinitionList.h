#ifdef DefineSSDParameter

// Both Building and Searching
DefineSSDParameter(m_buildSsdIndex, bool, true, "BuildSsdIndex")
DefineSSDParameter(m_vectorIDTranslate, string, string(""), "VectorIDTranslate")
DefineSSDParameter(m_headIndexFolder, string, string(""), "HeadIndexFolder")
DefineSSDParameter(m_queryFile, string, string(""), "QueryFile")
DefineSSDParameter(m_internalResultNum, int, 0, "InternalResultNum")
DefineSSDParameter(m_iNumberOfThreads, int, 16L, "NumberOfThreads")
DefineSSDParameter(m_headConfig, string, string(""), "HeadConfig") // must be in "Index" section

// Building
DefineSSDParameter(m_ssdIndex, string, string(""), "SsdIndex")
DefineSSDParameter(m_replicaCount, int, 4, "ReplicaCount")
DefineSSDParameter(m_postingPageLimit, int, 2, "PostingPageLimit")
DefineSSDParameter(m_outputEmptyReplicaID, bool, false, "OutputEmptyReplicaID")

// Searching
DefineSSDParameter(m_searchResult, string, string(""), "SearchResult")
DefineSSDParameter(m_extraFullGraphFile, string, string(""), "ExtraFullGraphFile")
DefineSSDParameter(m_extraGraphFile, string, string(""), "ExtraGraphFile")
DefineSSDParameter(m_extraGraphVectorSetFile, string, string(""), "ExtraGraphVectorSetFile")
DefineSSDParameter(m_truthFile, string, string(""), "TruthFile")
DefineSSDParameter(m_warmupFile, string, string(""), "WarmupFile")
DefineSSDParameter(m_logFile, string, string(""), "LogFile")
DefineSSDParameter(m_extraMaxCheck, string, string("1024"), "ExtraMaxCheck")
DefineSSDParameter(m_qpsLimit, int, 0, "QpsLimit")
DefineSSDParameter(m_resultNum, int, 32, "ResultNum")
DefineSSDParameter(m_queryCountLimit, int, (numeric_limits<int>::max)(), "QueryCountLimit")

#endif