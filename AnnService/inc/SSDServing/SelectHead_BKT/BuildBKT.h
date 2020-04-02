#pragma once

#include "inc/SSDServing/Common/stdafx.h"

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
			shared_ptr<COMMON::BKTree> BuildBKT(BasicVectorSet& p_vectorSet, const Options& opts) {
				HeadSimp<T> simp(p_vectorSet, opts.m_iDistCalcMethod);
				shared_ptr<COMMON::BKTree> bkt = make_shared<COMMON::BKTree>();
				bkt->m_iBKTKmeansK = opts.m_iBKTKmeansK;
				bkt->m_iBKTLeafSize = opts.m_iBKTLeafSize;
				bkt->m_iSamples = opts.m_iSamples;
				bkt->m_iTreeNumber = opts.m_iTreeNumber;
				fprintf(stdout, "Start invoking BuildTrees.\n");
				fprintf(stdout, "BKTKmeansK: %d, BKTLeafSize: %d, Samples: %d, TreeNumber: %d, ThreadNum: %d.\n", 
					bkt->m_iBKTKmeansK, bkt->m_iBKTLeafSize, bkt->m_iSamples, bkt->m_iTreeNumber, opts.m_iNumberOfThreads);
				auto start = std::chrono::system_clock::now();
				bkt->BuildTrees<T>(&simp, nullptr, nullptr, opts.m_iNumberOfThreads);
				auto end = std::chrono::system_clock::now();
				std::chrono::minutes elapsedMinutes = std::chrono::duration_cast<std::chrono::minutes>(end - start);
				fprintf(stdout, "End invoking BuildTrees.\n");
				fprintf(stdout, "Invoking BuildTrees used time: %d minutes (about %.2lf hours).\n", elapsedMinutes.count(), elapsedMinutes.count() / 60.0);

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
