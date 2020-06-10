#include "inc/SSDServing/SelectHead_BKT/BootSelectHead.h"
#include "inc/SSDServing/SelectHead_BKT/BuildBKT.h"
#include "inc/SSDServing/SelectHead_BKT/AnalyzeTree.h"
#include "inc/SSDServing/SelectHead_BKT/SelectHead.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"
#include "inc/SSDServing/IndexBuildManager/CommonDefines.h"

namespace SPTAG {
	namespace SSDServing {
		namespace SelectHead_BKT {
			ErrorCode Bootstrap(Options& opts) {

				VectorSearch::TimeUtils::StopW sw;

				fprintf(stdout, "Start loading vector file.\n");
				BasicVectorSet vectorSet(
					COMMON_OPTS.m_vectorPath.c_str(), 
					COMMON_OPTS.m_valueType, 
					COMMON_OPTS.m_dim, 
					COMMON_OPTS.m_vectorSize,
					COMMON_OPTS.m_vectorType,
					COMMON_OPTS.m_vectorDelimiter,
					COMMON_OPTS.m_distCalcMethod);
				fprintf(stdout, "Finish loading vector file.\n");

				fprintf(stdout, "Start generating BKT.\n");
				std::shared_ptr<COMMON::BKTree> bkt;
				switch (COMMON_OPTS.m_valueType)
				{
#define DefineVectorValueType(Name, Type) \
    case VectorValueType::Name: \
        bkt = BuildBKT<Type>(vectorSet, opts); \
		break;

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType

				default: break;
				}
				fprintf(stdout, "Finish generating BKT.\n");

				std::unordered_map<int, int> counter;

				if (opts.m_calcStd)
				{
					CalcLeafSize(0, bkt, counter);
				}

				if (opts.m_analyzeOnly)
				{
					fprintf(stdout, "Analyze Only.\n");

					std::vector<BKTNodeInfo> bktNodeInfos(bkt->size());

					// Always use the first tree
					DfsAnalyze(0, bkt, vectorSet, opts, 0, bktNodeInfos);

					fprintf(stdout, "Analyze Finish.\n");
				}
				else {
					if (SelectHead(vectorSet, bkt, opts, counter) != ErrorCode::Success)
					{
						return ErrorCode::Fail;
					}
				}

				double elapsedMinutes = sw.getElapsedMin();
				fprintf(stderr, "Total used time: %.2lf minutes (about %.2lf hours).\n", elapsedMinutes, elapsedMinutes / 60.0);

				return ErrorCode::Success;
			}
		}
	}
}