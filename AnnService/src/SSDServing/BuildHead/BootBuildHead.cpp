#include "inc/Core/Common.h"
#include "inc/Core/VectorIndex.h"

#include "inc/SSDServing/BuildHead/BootBuildHead.h"
#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"


namespace SPTAG {
	namespace SSDServing {
		namespace BuildHead {
			ErrorCode Bootstrap(Options& options, const SPTAG::Helper::IniReader::ParameterValueMap& params) {
                // These three params are mandatory.
                auto indexBuilder = SPTAG::VectorIndex::CreateInstance(COMMON_OPTS.m_indexAlgoType, COMMON_OPTS.m_valueType);
                indexBuilder->SetParameter("DistCalcMethod", SPTAG::Helper::Convert::ConvertToString(COMMON_OPTS.m_distCalcMethod));

                for (const auto& iter : params)
                {
                    indexBuilder->SetParameter(iter.first.c_str(), iter.second.c_str());
                }

                SPTAG::ErrorCode code;

                {
                    std::vector<std::string> files = SPTAG::Helper::StrUtils::SplitString(COMMON_OPTS.m_headVectorFile, ",");
                    std::ifstream inputStream(files[0], std::ifstream::binary);
                    if (!inputStream.is_open()) {
                        fprintf(stderr, "Failed to read input file.\n");
                        return ErrorCode::Fail;
                    }
                    SPTAG::SizeType row;
                    SPTAG::DimensionType col;
                    inputStream.read((char*)&row, sizeof(SPTAG::SizeType));
                    inputStream.read((char*)&col, sizeof(SPTAG::DimensionType));
                    std::uint64_t totalRecordVectorBytes = ((std::uint64_t)GetValueTypeSize(COMMON_OPTS.m_valueType)) * row * col;
                    SPTAG::ByteArray vectorSet = SPTAG::ByteArray::Alloc(totalRecordVectorBytes);
                    char* vecBuf = reinterpret_cast<char*>(vectorSet.Data());
                    inputStream.read(vecBuf, totalRecordVectorBytes);
                    inputStream.close();
                    std::shared_ptr<SPTAG::VectorSet> p_vectorSet(new SPTAG::BasicVectorSet(vectorSet, COMMON_OPTS.m_valueType, col, row));

                    std::shared_ptr<SPTAG::MetadataSet> p_metaSet = nullptr;
                    if (files.size() >= 3) {
                        p_metaSet.reset(new SPTAG::FileMetadataSet(files[1], files[2]));
                    }
                    code = indexBuilder->BuildIndex(p_vectorSet, p_metaSet);
                    indexBuilder->SaveIndex(COMMON_OPTS.m_headIndexFolder);
                }

                if (SPTAG::ErrorCode::Success != code)
                {
                    fprintf(stderr, "Failed to build index.\n");
                    return ErrorCode::Fail;
                }
				return ErrorCode::Success;
			}
		}
	}
}
