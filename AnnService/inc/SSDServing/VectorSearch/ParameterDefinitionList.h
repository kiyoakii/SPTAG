#ifdef DefineSSDParameter

// Both Building and Searching
DefineSSDParameter(m_buildSsdIndex, bool, false, "BuildSsdIndex")
DefineSSDParameter(m_vectorIDTranslate, std::string, std::string(""), "VectorIDTranslate")
DefineSSDParameter(m_headIndexFolder, std::string, std::string(""), "HeadIndexFolder")
DefineSSDParameter(m_internalResultNum, int, 0, "InternalResultNum")
DefineSSDParameter(m_iNumberOfThreads, int, 16L, "NumberOfThreads")
DefineSSDParameter(m_headConfig, std::string, std::string(""), "HeadConfig") // must be in "Index" section
DefineSSDParameter(m_ssdIndex, std::string, std::string(""), "SsdIndex")

// Building
DefineSSDParameter(m_replicaCount, int, 4, "ReplicaCount")
DefineSSDParameter(m_postingPageLimit, int, 2, "PostingPageLimit")
DefineSSDParameter(m_outputEmptyReplicaID, bool, false, "OutputEmptyReplicaID")

// Searching
DefineSSDParameter(m_searchResult, std::string, std::string(""), "SearchResult")
DefineSSDParameter(m_logFile, std::string, std::string(""), "LogFile")
DefineSSDParameter(m_qpsLimit, int, 0, "QpsLimit")
DefineSSDParameter(m_resultNum, int, 32, "ResultNum")
DefineSSDParameter(m_queryCountLimit, int, (std::numeric_limits<int>::max)(), "QueryCountLimit")

#endif