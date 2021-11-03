#include <climits>
#ifdef DefineSSDParameter

// Both Building and Searching
DefineSSDParameter(m_execute, bool, false, "isExecute")
DefineSSDParameter(m_buildSsdIndex, bool, false, "BuildSsdIndex")
DefineSSDParameter(m_internalResultNum, int, 64, "InternalResultNum")
DefineSSDParameter(m_iNumberOfThreads, int, 16, "NumberOfThreads")
DefineSSDParameter(m_headConfig, std::string, std::string(""), "HeadConfig") // must be in "Index" section

// Building
DefineSSDParameter(m_replicaCount, int, 8, "ReplicaCount")
DefineSSDParameter(m_postingPageLimit, int, 3, "PostingPageLimit")
DefineSSDParameter(m_outputEmptyReplicaID, bool, false, "OutputEmptyReplicaID")

// Searching
DefineSSDParameter(m_searchResult, std::string, std::string(""), "SearchResult")
DefineSSDParameter(m_logFile, std::string, std::string(""), "LogFile")
DefineSSDParameter(m_qpsLimit, int, 0, "QpsLimit")
DefineSSDParameter(m_resultNum, int, 5, "ResultNum")
DefineSSDParameter(m_maxCheck, int, 4096, "MaxCheck")
DefineSSDParameter(m_queryCountLimit, int, (std::numeric_limits<int>::max)(), "QueryCountLimit")
DefineSSDParameter(m_indexSize, double, 0.05, "IndexSize")
DefineSSDParameter(m_insertVectorsPath, std::string, std::string(""), "InsertVectorsPath")
DefineSSDParameter(m_headDistPostingnum, std::string, std::string(""), "HeadDistPostingnum")
DefineSSDParameter(m_randomDisabled, bool, false, "RandomDisabled")
DefineSSDParameter(m_SSDVectorDistPath, std::string, std::string(""), "SSDVectorDistPath")

DefineSSDParameter(m_truthFilePrefix, std::string, std::string(""), "TruthFilePrefix")
DefineSSDParameter(m_step, int, 0, "Step")
DefineSSDParameter(m_k, int, 2, "ClusterNum")
DefineSSDParameter(m_searchVectorLimit, int, INT_MAX, "SearchVectorLimit")
DefineSSDParameter(m_storeSearchDetail, bool, false, "StoreSearchDetail")
#endif