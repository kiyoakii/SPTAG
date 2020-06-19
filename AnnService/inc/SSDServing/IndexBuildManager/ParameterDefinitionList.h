// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.
#ifdef DefineBasicParameter

// DefineBasicParameter(VarName, VarType, DefaultValue, RepresentStr)
DefineBasicParameter(m_valueType, SPTAG::VectorValueType, SPTAG::VectorValueType::Undefined, "ValueType")
DefineBasicParameter(m_distCalcMethod, SPTAG::DistCalcMethod, SPTAG::DistCalcMethod::Undefined, "DistCalcMethod")
DefineBasicParameter(m_dim, SPTAG::DimensionType, -1, "Dim")
DefineBasicParameter(m_vectorPath, std::string, std::string(""), "VectorPath")
DefineBasicParameter(m_vectorType, SPTAG::VectorFileType, SPTAG::VectorFileType::Undefined, "VectorType")
DefineBasicParameter(m_vectorSize, SPTAG::SizeType, -1, "VectorSize")
DefineBasicParameter(m_vectorDelimiter, std::string, std::string("|"), "VectorDelimiter")
DefineBasicParameter(m_queryPath, std::string, std::string(""), "QueryPath")
DefineBasicParameter(m_queryType, SPTAG::VectorFileType, SPTAG::VectorFileType::Undefined, "QueryType")
DefineBasicParameter(m_querySize, SPTAG::SizeType, -1, "QuerySize")
DefineBasicParameter(m_queryDelimiter, std::string, std::string("|"), "QueryDelimiter")
DefineBasicParameter(m_warmupPath, std::string, std::string(""), "WarmupPath")
DefineBasicParameter(m_warmupType, SPTAG::VectorFileType, SPTAG::VectorFileType::Undefined, "WarmupType")
DefineBasicParameter(m_warmupSize, SPTAG::SizeType, -1, "WarmupSize")
DefineBasicParameter(m_warmupDelimiter, std::string, std::string("|"), "WarmupDelimiter")
DefineBasicParameter(m_truthPath, std::string, std::string(""), "TruthPath")
DefineBasicParameter(m_truthType, SPTAG::TruthFileType, SPTAG::TruthFileType::Undefined, "TruthType")

#endif
