// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
#ifdef DefineBuildHeadParameter

// DefineBuildHeadParameter(VarName, VarType, DefaultValue, RepresentStr)
DefineBuildHeadParameter(m_inputFiles, std::string, std::string("vectors.bin"), "HeadVectorFile")
DefineBuildHeadParameter(m_outputFolder, std::string, std::string("HeadVectors.bin"), "HeadIndex")
DefineBuildHeadParameter(m_indexAlgoType, SPTAG::IndexAlgoType, SPTAG::IndexAlgoType::BKT, "IndexAlgoType")
DefineBuildHeadParameter(m_builderConfigFile, std::string, std::string("builder.ini"), "BuilderConfigFile")
DefineBuildHeadParameter(m_threadNum, std::uint32_t, omp_get_num_threads(), "NumberOfThreads")

#endif
