#pragma once

#ifdef DefineSSDParameter

DefineSSDParameter(m_buildSsdIndex, bool, true, "BuildSsdIndex")
DefineSSDParameter(m_extraMaxCheck, string, "1024", "ExtraMaxCheck")
DefineSSDParameter(m_parallelLoadPercentage, string, "0", "ParallelLoadPercentage")
DefineSSDParameter(m_iNumberOfThreads, int, 16L, "NumberOfThreads")

#endif