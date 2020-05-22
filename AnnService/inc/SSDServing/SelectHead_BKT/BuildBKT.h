#pragma once
#include "inc/Core/Common/BKTree.h"
#include "inc/Core/Common/DistanceUtils.h"
#include "inc/SSDServing/VectorSearch/TimeUtils.h"

namespace SPTAG {
	namespace SSDServing {
		namespace SelectHead_BKT {
			template<typename T>
			class HeadSimp : public COMMON::Simp {
			private:
				BasicVectorSet& m_vectorSet;
				DistCalcMethod m_iDistCalcMethod;
				float(*m_fComputeDistance)(const T* pX, const T* pY, DimensionType length);
			public:
				HeadSimp(BasicVectorSet& p_vectorSet, DistCalcMethod p_distCalcMethod) :
					m_vectorSet(p_vectorSet),
					m_iDistCalcMethod(p_distCalcMethod),
					m_fComputeDistance(COMMON::DistanceCalcSelector<T>(m_iDistCalcMethod))
				{
				}

				~HeadSimp() {}

				SizeType GetNumSamples() {
					return m_vectorSet.Count();
				}

				DimensionType GetFeatureDim() {
					return m_vectorSet.Dimension();
				}

				float ComputeDistance(const void* pX, const void* pY) {
					return m_fComputeDistance((const T*)pX, (const T*)pY, m_vectorSet.Dimension());
				}

				DistCalcMethod GetDistCalcMethod() {
					return m_iDistCalcMethod;
				}

				const void* GetSample(const SizeType idx) {
					return m_vectorSet.GetVector(idx);
				}

			};

			template<typename T>
			std::shared_ptr<COMMON::BKTree> BuildBKT(BasicVectorSet& p_vectorSet, const Options& opts) {
				HeadSimp<T> simp(p_vectorSet, opts.m_iDistCalcMethod);
				std::shared_ptr<COMMON::BKTree> bkt = std::make_shared<COMMON::BKTree>();
				bkt->m_iBKTKmeansK = opts.m_iBKTKmeansK;
				bkt->m_iBKTLeafSize = opts.m_iBKTLeafSize;
				bkt->m_iSamples = opts.m_iSamples;
				bkt->m_iTreeNumber = opts.m_iTreeNumber;
				fprintf(stdout, "Start invoking BuildTrees.\n");
				fprintf(stdout, "BKTKmeansK: %d, BKTLeafSize: %d, Samples: %d, TreeNumber: %d, ThreadNum: %d.\n", 
					bkt->m_iBKTKmeansK, bkt->m_iBKTLeafSize, bkt->m_iSamples, bkt->m_iTreeNumber, opts.m_iNumberOfThreads);
				VectorSearch::TimeUtils::StopW sw;
				bkt->BuildTrees<T>(&simp, nullptr, nullptr, opts.m_iNumberOfThreads);
				double elapsedMinutes = sw.getElapsedMin();

				fprintf(stdout, "End invoking BuildTrees.\n");
				fprintf(stdout, "Invoking BuildTrees used time: %.2lf minutes (about %.2lf hours).\n", elapsedMinutes, elapsedMinutes / 60.0);

				std::stringstream bktFileNameBuilder;
				bktFileNameBuilder << opts.m_vectorFile << ".bkt."
					<< opts.m_iBKTKmeansK << "_"
					<< opts.m_iBKTLeafSize << "_"
					<< opts.m_iTreeNumber << "_"
					<< opts.m_iSamples << "_"
					<< static_cast<int>(opts.m_iDistCalcMethod) << ".bin";

				std::string bktFileName = bktFileNameBuilder.str();
				if (opts.m_saveBKT) {
					bkt->SaveTrees(bktFileName);
				}

				return bkt;
			}
		}
	}
}
