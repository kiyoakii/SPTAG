#include "inc/Test.h"
#include <string>
#include <fstream>
#include <atlstr.h>
#include <random>
#include <type_traits>
#include <functional>
#include <algorithm>

#include <boost/filesystem.hpp>

#include "inc/Core/Common.h"
#include "inc/Helper/StringConvert.h"
#include "inc/SSDServing/IndexBuildManager/main.h"
#include "inc/Core/Common/DistanceUtils.h"
#include "inc/Core/Common/CommonUtils.h"
#include <ppl.h>

using namespace std;

template<typename T>
void GenerateVectors(string fileName, SPTAG::SizeType rows, SPTAG::DimensionType dims, SPTAG::VectorFileType fileType) {
	if (boost::filesystem::exists(fileName))
	{
		fprintf(stdout, "%s was generated. Skip generation.", fileName.c_str());
		return;
	}
	
	ofstream of(fileName, ofstream::binary);
	if (!of.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", fileName.c_str());
		BOOST_CHECK(false);
		return;
	}

	uniform_real_distribution<float> ud(0, 126);
	mt19937 mt(543);
	vector<T> tmp(dims);

	if (fileType == SPTAG::VectorFileType::DEFAULT)
	{
		of.write(reinterpret_cast<char*>(&rows), 4);
		of.write(reinterpret_cast<char*>(&dims), 4);

		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < dims; j++)
			{
				float smt = ud(mt);
				tmp[j] = static_cast<T>(smt);
			}

			SPTAG::COMMON::Utils::Normalize(tmp.data(), dims, SPTAG::COMMON::Utils::GetBase<T>());
			of.write(reinterpret_cast<char*>(tmp.data()), dims * sizeof(T));
		}
	}
	else if (fileType == SPTAG::VectorFileType::XVEC)
	{
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < dims; j++)
			{
				float smt = ud(mt);
				tmp[j] = static_cast<T>(smt);
			}

			SPTAG::COMMON::Utils::Normalize(tmp.data(), dims, SPTAG::COMMON::Utils::GetBase<T>());
			of.write(reinterpret_cast<char*>(&dims), 4);
			of.write(reinterpret_cast<char*>(tmp.data()), dims * sizeof(T));
		}

	}

}

void GenVec(string vectorsName, SPTAG::VectorValueType vecType, SPTAG::VectorFileType vecFileType, SPTAG::SizeType rows = 1000, SPTAG::DimensionType dims = 100) {
	switch (vecType)
	{
#define DefineVectorValueType(Name, Type) \
case SPTAG::VectorValueType::Name: \
GenerateVectors<Type>(vectorsName, rows, dims, vecFileType); \
break; \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
	default:
		break;
	}
}

template<typename T>
void FillVectors(const string binName, vector<vector<T>>& bins, SPTAG::VectorFileType vft) {
	ifstream f(binName, ifstream::binary);
	if (!f.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", binName.c_str());
		BOOST_CHECK(false);
		return;
	}

	if (vft == SPTAG::VectorFileType::DEFAULT)
	{
		SPTAG::SizeType rowNum;
		SPTAG::DimensionType dimNum;
		f.read(reinterpret_cast<char*>(&rowNum), sizeof(SPTAG::SizeType));
		f.read(reinterpret_cast<char*>(&dimNum), sizeof(SPTAG::DimensionType));
		vector<T> cur(dimNum);
		for (size_t i = 0; i < rowNum; i++)
		{
			f.read(reinterpret_cast<char*>(cur.data()), dimNum * sizeof(T));
			bins.push_back(cur);
		}
		
	}
	else if(vft == SPTAG::VectorFileType::XVEC)
	{
		vector<T> cur;
		SPTAG::DimensionType dimNum;
		do {
			f.read(reinterpret_cast<char*>(&dimNum), 4);
			if (f.eof()) break;
			if (!f.good())
			{
				fprintf(stderr, "ERROR: file %s isn't good.\n", binName.c_str());
				exit(1);
			}
			cur.resize(dimNum);
			f.read(reinterpret_cast<char*>(cur.data()), dimNum * sizeof(T));
			bins.push_back(cur);
		} while (true);

	}
}

struct Neighbor
{
	SPTAG::SizeType key;
	float dist;

	Neighbor(SPTAG::SizeType k, float d) : key(k), dist(d) {}

	bool operator < (const Neighbor& another) const
	{
		return this->dist == another.dist ? this->key < another.key : this->dist < another.dist;
	}
};

void writeTruthFile(const string truthFile, size_t queryNumber, const int K, vector<vector<SPTAG::SizeType>>& truthset, SPTAG::TruthFileType TFT) {
	
	if (TFT == SPTAG::TruthFileType::TXT)
	{
		ofstream of(truthFile);
		for (size_t i = 0; i < queryNumber; i++)
		{
			for (size_t k = 0; k < K; k++)
			{
				of << truthset[i][k];
				if (k != K - 1)
				{
					of << " ";
				}
			}
			of << endl;
		}
	}
	else if (TFT == SPTAG::TruthFileType::XVEC)
	{
		ofstream of(truthFile, ios_base::binary);
		for (size_t i = 0; i < queryNumber; i++)
		{
			of.write(reinterpret_cast<const char*>(&K), 4);
			of.write(reinterpret_cast<char*>(truthset[i].data()), K * 4);
		}
	}
}

template<typename T>
void GenerateTruth(const string queryFile, const string vectorFile, const string truthFile,
	const SPTAG::DistCalcMethod distMethod, const int K, SPTAG::VectorFileType p_vecFileType, const SPTAG::TruthFileType p_truthFileType) {
	if (boost::filesystem::exists(truthFile))
	{
		fprintf(stdout, "truthFile: %s was generated. Skip generation.", truthFile.c_str());
		return;
	}
	vector<vector<T>> querys;
	vector<vector<T>> vectors;
	FillVectors<T>(queryFile, querys, p_vecFileType);
	FillVectors<T>(vectorFile, vectors, p_vecFileType);
	vector<vector<SPTAG::SizeType>> truthset(querys.size(), vector<SPTAG::SizeType>(K, 0));
	std::atomic_uint32_t processed = 0;

	concurrency::parallel_for(0, 16, [&](int tid)
		{
			// LARGE_INTEGER timePoint;
			for (uint32_t i = processed.fetch_add(1); i < querys.size(); i = processed.fetch_add(1))
			{
				vector<T> curQuery = querys[i];
				vector<Neighbor> neighbours;
				bool isFirst = true;
				for (size_t j = 0; j < vectors.size(); j++)
				{
					vector<T> curVector = vectors[j];
					if (curQuery.size() != curVector.size())
					{
						fprintf(stderr, "query and vector have different dimensions.");
						BOOST_CHECK(false);
						return;
					}

					float dist;

					if (distMethod == SPTAG::DistCalcMethod::L2)
					{
						dist = SPTAG::COMMON::DistanceUtils::ComputeL2Distance(curQuery.data(), curVector.data(), curQuery.size());
					}
					else {
						dist = SPTAG::COMMON::DistanceUtils::ComputeCosineDistance(curQuery.data(), curVector.data(), curQuery.size());
					}

					Neighbor nei(j, dist);
					neighbours.push_back(nei);
					if (neighbours.size() == K && isFirst)
					{
						make_heap(neighbours.begin(), neighbours.end());
						isFirst = false;
					}
					if (neighbours.size() > K)
					{
						push_heap(neighbours.begin(), neighbours.end());
						pop_heap(neighbours.begin(), neighbours.end());
						neighbours.pop_back();
					}
				}

				if (K != neighbours.size())
				{
					fprintf(stderr, "K is too big.\n");
					BOOST_CHECK(false);
					return;
				}

				std::sort(neighbours.begin(), neighbours.end());

				for (size_t k = 0; k < K; k++)
				{
					truthset[i][k] = neighbours[k].key;
				}

			}
		});
	
	writeTruthFile(truthFile, querys.size(), K, truthset, p_truthFileType);
}

void GenerateTruth(const string queryFile, const string vectorFile, const string truthFile,
	const SPTAG::DistCalcMethod distMethod, const int K, SPTAG::VectorValueType vvt, SPTAG::VectorFileType p_vecFileType, SPTAG::TruthFileType p_truthFileType) {
#define DefineVectorValueType(Name, Type) \
	if (vvt == SPTAG::VectorValueType::Name) { \
		GenerateTruth<Type>(queryFile, vectorFile, truthFile, distMethod, K, p_vecFileType, p_truthFileType); \
	} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
}

void TestHead(string vectorsName, string configName, string OutputIDFile, string OutputVectorFile, 
	SPTAG::VectorValueType vecType, SPTAG::DistCalcMethod distMethod, 
	SPTAG::VectorFileType p_vectorFileType, SPTAG::SizeType p_iVectorNumber, SPTAG::DimensionType p_iDimension,
	int p_BKTKmeansK, int p_BKTLeafSize, int p_SamplesNumber,
	int p_SelectThreshold, int p_SplitFactor, int p_SplitThreshold) {

	ofstream config(configName);
	if (!config.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", configName.c_str());
		BOOST_CHECK(false);
		return;
	}
	config << "[SelectHead]" << endl;
	config << "VectorFilePath=" << vectorsName << endl;
	config << "VectorValueType=" << SPTAG::Helper::Convert::ConvertToString(vecType) << endl;
	config << "DistCalcMethod=" << SPTAG::Helper::Convert::ConvertToString(distMethod) << endl;

	config << "VectorFileType=" << SPTAG::Helper::Convert::ConvertToString(p_vectorFileType) << endl;
	config << "VectorNumber=" << p_iVectorNumber << endl;
	config << "Dimension=" << p_iDimension << endl;

	config << "TreeNumber=" << "1" << endl;
	config << "BKTKmeansK=" << p_BKTKmeansK << endl;
	config << "BKTLeafSize=" << p_BKTLeafSize << endl;
	config << "SamplesNumber=" << p_SamplesNumber << endl;
	config << "NumberOfThreads=" << "1" << endl;
	config << "SaveBKT=" << "true" <<endl;

	config << "AnalyzeOnly=" << "false" << endl;
	config << "CalcStd=" << "true" << endl;
	config << "SelectDynamically=" << "true" << endl;
	config << "NoOutput=" << "false" << endl;

	config << "SelectThreshold=" << p_SelectThreshold << endl;
	config << "SplitFactor=" << p_SplitFactor << endl;
	config << "SplitThreshold=" << p_SplitThreshold << endl;
	config << "Ratio=" << "0.2" << endl;
	config << "RecursiveCheckSmallCluster=" << "true" << endl;
	config << "PrintSizeCount=" << "true" << endl;

	config << "OutputIDFile=" << OutputIDFile << endl;
	config << "OutputVectorFile=" << OutputVectorFile << endl;

	config.close();

	char* arg1 = new char[100];
	char* arg2 = new char[100];
	strcpy_s(arg1, 100, "SSDServing.exe");
	strcpy_s(arg2, 100, configName.c_str());
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
	delete[] arg1;
	delete[] arg2;
}

void TestBuildHead(
	string configName, string p_headVectorFile, string p_headIndexFile,
	SPTAG::IndexAlgoType p_indexAlgoType, SPTAG::DistCalcMethod p_distCalcMethod,
	string p_builderFile, SPTAG::VectorValueType p_inputValueType, 
	int p_BKTKmeansK, int p_BKTLeafSize, int p_Samples,
	uint32_t p_threadNum) {

	ofstream config(configName);
	if (!config.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", configName.c_str());
		BOOST_CHECK(false);
		return;
	}

	{
		if (p_builderFile.empty())
		{
			fprintf(stderr, "no builder file for head index build.\n");
			BOOST_CHECK(false);
			return;
		}

		if (boost::filesystem::exists(p_builderFile))
		{
			fprintf(stdout, "%s was generated. Skip generation.", p_builderFile.c_str());
		}
		else {
			ofstream bf(p_builderFile);
			if (!bf.is_open())
			{
				fprintf(stderr, "%s can't be opened. ", p_builderFile.c_str());
				BOOST_CHECK(false);
				return;
			}
			bf << "[Index]" << endl;
			bf << "DistCalcMethod=" << SPTAG::Helper::Convert::ConvertToString(p_distCalcMethod) << endl;
			bf << "BKTKmeansK=" << p_BKTKmeansK << endl;
			bf << "BKTLeafSize=" << p_BKTLeafSize << endl;
			bf << "Samples=" << p_Samples << endl;
			bf.close();
		}
	}

	config << "[BuildHead]" << endl;
	config << "HeadVectorFile=" << p_headVectorFile << endl;
	config << "HeadIndex=" << p_headIndexFile << endl;
	config << "IndexAlgoType=" << SPTAG::Helper::Convert::ConvertToString(p_indexAlgoType) << endl;
	config << "BuilderConfigFile=" << p_builderFile << endl;
	config << "VectorValueType=" << SPTAG::Helper::Convert::ConvertToString(p_inputValueType) << endl;
	config << "ThreadNum=" << p_threadNum << endl;

	config.close();

	char* arg1 = new char[100];
	char* arg2 = new char[100];
	strcpy_s(arg1, 100, "SSDServing.exe");
	strcpy_s(arg2, 100, configName.c_str());
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
	delete[] arg1;
	delete[] arg2;
}

void TestBuildSSDIndex(string configName,
	string p_vectorIDTranslate,
	string p_headIndexFolder,
	string p_queryFile,
	string p_ssdIndex,
	string p_headConfig,

	SPTAG::VectorFileType p_queryFileType,
	SPTAG::SizeType p_iQueryNumber,
	SPTAG::DimensionType p_iQueryDimension,

	int p_internalResultNum,
	int p_numberOfThreads,
	int p_replicaCount,
	int p_postingPageLimit,
	bool p_outputEmptyReplicaID
){
	ofstream config(configName);
	if (!config.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", configName.c_str());
		BOOST_CHECK(false);
		return;
	}

	if (!p_headConfig.empty())
	{
		if (boost::filesystem::exists(p_headConfig))
		{
			fprintf(stdout, "%s was generated. Skip generation.", p_headConfig.c_str());
		}
		else {
			ofstream bf(p_headConfig);
			if (!bf.is_open())
			{
				fprintf(stderr, "%s can't be opened. ", p_headConfig.c_str());
				BOOST_CHECK(false);
				return;
			}
			bf << "[Index]" << endl;
			bf.close();
		}
	}

	config << "[BuildSSDIndex]" << endl;
	config << "BuildSsdIndex=" << "true" << endl;
	config << "VectorIDTranslate=" << p_vectorIDTranslate << endl;
	config << "HeadIndexFolder=" << p_headIndexFolder << endl;
	config << "QueryFile=" << p_queryFile << endl;
	config << "InternalResultNum=" << p_internalResultNum << endl;
	config << "NumberOfThreads=" << p_numberOfThreads << endl;
	config << "HeadConfig=" << p_headConfig << endl;

	config << "QueryFileType=" << SPTAG::Helper::Convert::ConvertToString(p_queryFileType) << endl;
	config << "QueryNumber=" << p_iQueryNumber << endl;
	config << "QueryDimension=" << p_iQueryDimension << endl;
	
	config << "SsdIndex=" << p_ssdIndex << endl;
	config << "ReplicaCount=" << p_replicaCount << endl;
	config << "PostingPageLimit=" << p_postingPageLimit << endl;
	config << "OutputEmptyReplicaID=" << p_outputEmptyReplicaID << endl;

	config.close();

	char* arg1 = new char[100];
	char* arg2 = new char[100];
	strcpy_s(arg1, 100, "SSDServing.exe");
	strcpy_s(arg2, 100, configName.c_str());
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
	delete[] arg1;
	delete[] arg2;
}

void TestSearchSSDIndex(string configName,
	string p_vectorIDTranslate, string p_headIndexFolder, string p_extraFullGraphFile, string p_headConfig, string p_queryFile,
	int p_internalResultNum, int p_resultNum, int p_numberOfThreads,

	string p_searchResult,
	string p_extraGraphFile,
	string p_extraGraphVectorSetFile,
	string p_truthFile,
	string p_warmupFile,
	string p_logFile,
	string p_ExtraMaxCheck,
	int p_qpsLimit,
	int p_queryCountLimit,

	SPTAG::VectorFileType p_queryFileType, SPTAG::SizeType p_iQueryNumber, SPTAG::DimensionType p_iQueryDimension,
	SPTAG::VectorFileType p_warmupFileType, SPTAG::SizeType p_iWarmupNumber, SPTAG::DimensionType p_iWarmupDimension,
	SPTAG::TruthFileType p_truthFileType, SPTAG::SizeType p_iTruthNumber
) {
	ofstream config(configName);
	if (!config.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", configName.c_str());
		BOOST_CHECK(false);
		return;
	}

	if (!p_headConfig.empty())
	{
		if (boost::filesystem::exists(p_headConfig))
		{
			fprintf(stdout, "%s was generated. Skip generation.", p_headConfig.c_str());
		}
		else {
			ofstream bf(p_headConfig);
			if (!bf.is_open())
			{
				fprintf(stderr, "%s can't be opened. ", p_headConfig.c_str());
				BOOST_CHECK(false);
				return;
			}
			bf << "[Index]" << endl;
			bf.close();
		}
	}

	config << "[SearchSSDIndex]" << endl;
	config << "BuildSsdIndex=" << "false" << endl;
	config << "VectorIDTranslate=" << p_vectorIDTranslate << endl;
	config << "HeadIndexFolder=" << p_headIndexFolder << endl;
	config << "QueryFile=" << p_queryFile << endl;
	config << "InternalResultNum=" << p_internalResultNum << endl;
	config << "NumberOfThreads=" << p_numberOfThreads << endl;

	config << "SearchResult=" << p_searchResult << endl;
	config << "ExtraFullGraphFile=" << p_extraFullGraphFile << endl;
	config << "HeadConfig=" << p_headConfig << endl;
	config << "ExtraGraphFile=" << p_extraGraphFile << endl;
	config << "ExtraGraphVectorSetFile=" << p_extraGraphVectorSetFile << endl;
	config << "TruthFile=" << p_truthFile << endl;
	config << "WarmupFile=" << p_warmupFile << endl;
	config << "LogFile=" << p_logFile << endl;
	config << "ExtraMaxCheck=" << p_ExtraMaxCheck << endl;
	config << "QpsLimit=" << p_qpsLimit << endl;
	config << "ResultNum=" << p_resultNum << endl;
	config << "QueryCountLimit=" << p_queryCountLimit << endl;


	config << "QueryFileType=" << SPTAG::Helper::Convert::ConvertToString(p_queryFileType) << endl;
	config << "QueryNumber=" << p_iQueryNumber << endl;
	config << "QueryDimension=" << p_iQueryDimension << endl;
	config << "WarmupFileType=" << SPTAG::Helper::Convert::ConvertToString(p_warmupFileType) << endl;
	config << "WarmupNumber=" << p_iWarmupNumber << endl;
	config << "WarmupDimension=" << p_iWarmupDimension << endl;
	config << "TruthFileType=" << SPTAG::Helper::Convert::ConvertToString(p_truthFileType) << endl;
	config << "TruthNumber=" << p_iTruthNumber << endl;

	config.close();

	char* arg1 = new char[100];
	char* arg2 = new char[100];
	strcpy_s(arg1, 100, "SSDServing.exe");
	strcpy_s(arg2, 100, configName.c_str());
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
	delete[] arg1;
	delete[] arg2;
}

BOOST_AUTO_TEST_SUITE(SSDServingTest)

#define SSDTEST_DIRECTORY_NAME "sddtest"
#define VECTOR_NUM 1000
#define QUERY_NUM 10
#define VECTOR_DIM 100
#define SSDTEST_DIRECTORY SSDTEST_DIRECTORY_NAME "\\"
#define VECTORS(VT, FT) SSDTEST_DIRECTORY "vectors_"#VT"_"#FT".bin"
#define QUERIES(VT, FT) SSDTEST_DIRECTORY "vectors_"#VT"_"#FT".query"
#define TRUTHSET(VT, DM, FT, TFT) SSDTEST_DIRECTORY "vectors_"#VT"_"#DM"_"#FT"_"#TFT".truth"
#define HEAD_IDS(VT, DM, FT) SSDTEST_DIRECTORY "head_ids_"#VT"_"#DM"_"#FT".bin"
#define HEAD_VECTORS(VT, DM, FT) SSDTEST_DIRECTORY "head_vectors_"#VT"_"#DM"_"#FT".bin"
#define HEAD_INDEX(VT, DM, ALGO, FT) SSDTEST_DIRECTORY "head_"#VT"_"#DM"_"#ALGO"_"#FT".head_index"
#define SSD_INDEX(VT, DM, ALGO, FT) SSDTEST_DIRECTORY "ssd_"#VT"_"#DM"_"#ALGO"_"#FT".ssd_index"

#define SELECT_HEAD_CONFIG(VT, DM, FT) SSDTEST_DIRECTORY "test_head_"#VT"_"#DM"_"#FT".ini"
#define BUILD_HEAD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_head_"#VT"_"#DM"_"#ALGO".ini"
#define BUILD_HEAD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_head_"#VT"_"#DM"_"#ALGO".builder.ini"
#define BUILD_SSD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_ssd"#VT"_"#DM"_"#ALGO".ini"
#define BUILD_SSD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_ssd"#VT"_"#DM"_"#ALGO".builder.ini"
#define SEARCH_SSD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO".ini"
#define SEARCH_SSD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO".builder.ini"
#define SEARCH_SSD_RESULT(VT, DM, ALGO, FT, TFT) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO"_"#FT"_"#TFT".result"

#define GVQ(VT, FT) \
BOOST_AUTO_TEST_CASE(GenerateVectorsQueries##VT##FT) { \
boost::filesystem::create_directory(SSDTEST_DIRECTORY_NAME); \
GenVec(VECTORS(VT, FT), SPTAG::VectorValueType::VT, SPTAG::VectorFileType::FT, VECTOR_NUM, VECTOR_DIM); \
GenVec(QUERIES(VT, FT), SPTAG::VectorValueType::VT, SPTAG::VectorFileType::FT, QUERY_NUM, VECTOR_DIM); \
} \

GVQ(Float, DEFAULT)
GVQ(Int16, DEFAULT)
GVQ(UInt8, DEFAULT)
GVQ(Int8, DEFAULT)

GVQ(Float, XVEC)
GVQ(Int16, XVEC)
GVQ(UInt8, XVEC)
GVQ(Int8, XVEC)
#undef GVQ

#define GTR(VT, DM, FT, TFT) \
BOOST_AUTO_TEST_CASE(GenerateTruth##VT##DM##FT##TFT) { \
GenerateTruth(QUERIES(VT, FT), VECTORS(VT, FT), TRUTHSET(VT, DM, FT, TFT), \
	SPTAG::DistCalcMethod::DM, 128, SPTAG::VectorValueType::VT, SPTAG::VectorFileType::FT, SPTAG::TruthFileType::TFT); \
} \

GTR(Float, L2, DEFAULT, TXT)
GTR(Float, Cosine, DEFAULT, TXT)
GTR(Int16, L2, DEFAULT, TXT)
GTR(Int16, Cosine, DEFAULT, TXT)
GTR(UInt8, L2, DEFAULT, TXT)
GTR(UInt8, Cosine, DEFAULT, TXT)
GTR(Int8, L2, DEFAULT, TXT)
GTR(Int8, Cosine, DEFAULT, TXT)

GTR(Float, L2, XVEC, XVEC)
GTR(Float, Cosine, XVEC, XVEC)
GTR(Int16, L2, XVEC, XVEC)
GTR(Int16, Cosine, XVEC, XVEC)
GTR(UInt8, L2, XVEC, XVEC)
GTR(UInt8, Cosine, XVEC, XVEC)
GTR(Int8, L2, XVEC, XVEC)
GTR(Int8, Cosine, XVEC, XVEC)
#undef GTR

#define WTEV(VT, DM, FT, \
	S_BKTKmeansK, S_BKTLeafSize, S_SamplesNumber, \
S_SelectThreshold, S_SplitFactor, S_SplitThreshold) \
BOOST_AUTO_TEST_CASE(TestHead##VT##DM##FT) \
{ \
	string vectorsName = VECTORS(VT, FT); \
	string configName = SELECT_HEAD_CONFIG(VT, DM, FT); \
	string OutputIDFile = HEAD_IDS(VT, DM, FT); \
	string OutputVectorFile = HEAD_VECTORS(VT, DM, FT); \
\
	TestHead(vectorsName, configName, OutputIDFile, OutputVectorFile, \
		SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, \
		SPTAG::VectorFileType::FT, VECTOR_NUM, VECTOR_DIM, \
		S_BKTKmeansK, S_BKTLeafSize, S_SamplesNumber, \
		S_SelectThreshold, S_SplitFactor, S_SplitThreshold \
	); \
} \

WTEV(Float, L2, DEFAULT,
	3, 6, 100, 
	12, 9, 18)
WTEV(Float, Cosine, DEFAULT,
	3, 6, 100, 
	12, 9, 18)
WTEV(Int16, L2, DEFAULT,
	3, 6, 100,
	12, 9, 18)
WTEV(Int16, Cosine, DEFAULT,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, L2, DEFAULT,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, Cosine, DEFAULT,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, L2, DEFAULT,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, Cosine, DEFAULT,
	3, 6, 100,
	12, 9, 18)


WTEV(Float, L2, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(Float, Cosine, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(Int16, L2, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(Int16, Cosine, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, L2, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, Cosine, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, L2, XVEC,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, Cosine, XVEC,
	3, 6, 100,
	12, 9, 18)

#undef WTEV

#define BDHD(VT, DM, ALGO, FT, \
H_BKTKmeansK, H_BKTLeafSize, H_Samples) \
BOOST_AUTO_TEST_CASE(TestBuildHead##VT##DM##ALGO##FT) { \
string configName = BUILD_HEAD_CONFIG(VT, DM, ALGO); \
string builderFile = BUILD_HEAD_BUILDER_CONFIG(VT, DM, ALGO); \
TestBuildHead( \
	configName, \
	HEAD_VECTORS(VT, DM, FT), \
	HEAD_INDEX(VT, DM, ALGO, FT), \
	SPTAG::IndexAlgoType::ALGO, \
	SPTAG::DistCalcMethod::DM, \
	builderFile, \
	SPTAG::VectorValueType::VT, \
	H_BKTKmeansK, H_BKTLeafSize, H_Samples, \
	3 \
);} \

BDHD(Float, L2, BKT, DEFAULT,
	2, 4, 20)
BDHD(Float, L2, KDT, DEFAULT,
	2, 4, 20)
BDHD(Float, Cosine, BKT, DEFAULT,
	2, 4, 20)
BDHD(Float, Cosine, KDT, DEFAULT,
	2, 4, 20)

BDHD(Int8, L2, BKT, DEFAULT,
	2, 4, 20)
BDHD(Int8, L2, KDT, DEFAULT,
	2, 4, 20)
BDHD(Int8, Cosine, BKT, DEFAULT,
	2, 4, 20)
BDHD(Int8, Cosine, KDT, DEFAULT,
	2, 4, 20)

BDHD(UInt8, L2, BKT, DEFAULT,
	2, 4, 20)
BDHD(UInt8, L2, KDT, DEFAULT,
	2, 4, 20)
BDHD(UInt8, Cosine, BKT, DEFAULT,
	2, 4, 20)
BDHD(UInt8, Cosine, KDT, DEFAULT,
	2, 4, 20)

BDHD(Int16, L2, BKT, DEFAULT,
	2, 4, 20)
BDHD(Int16, L2, KDT, DEFAULT,
	2, 4, 20)
BDHD(Int16, Cosine, BKT, DEFAULT,
	2, 4, 20)
BDHD(Int16, Cosine, KDT, DEFAULT,
	2, 4, 20)


//XVEC
BDHD(Float, L2, BKT, XVEC,
	2, 4, 20)
BDHD(Float, L2, KDT, XVEC,
	2, 4, 20)
BDHD(Float, Cosine, BKT, XVEC,
	2, 4, 20)
BDHD(Float, Cosine, KDT, XVEC,
	2, 4, 20)

BDHD(Int8, L2, BKT, XVEC,
	2, 4, 20)
BDHD(Int8, L2, KDT, XVEC,
	2, 4, 20)
BDHD(Int8, Cosine, BKT, XVEC,
	2, 4, 20)
BDHD(Int8, Cosine, KDT, XVEC,
	2, 4, 20)

BDHD(UInt8, L2, BKT, XVEC,
	2, 4, 20)
BDHD(UInt8, L2, KDT, XVEC,
	2, 4, 20)
BDHD(UInt8, Cosine, BKT, XVEC,
	2, 4, 20)
BDHD(UInt8, Cosine, KDT, XVEC,
	2, 4, 20)

BDHD(Int16, L2, BKT, XVEC,
	2, 4, 20)
BDHD(Int16, L2, KDT, XVEC,
	2, 4, 20)
BDHD(Int16, Cosine, BKT, XVEC,
	2, 4, 20)
BDHD(Int16, Cosine, KDT, XVEC,
	2, 4, 20)

#undef BDHD

#define BDSSD(VT, DM, ALGO, FT, \
	BS_internalResultNum, BS_numberOfThreads, BS_replicaCount, BS_postingPageLimit) \
BOOST_AUTO_TEST_CASE(TestBuildSSDIndex##VT##DM##ALGO##FT) { \
string configName = BUILD_SSD_CONFIG(VT, DM, ALGO); \
TestBuildSSDIndex(\
	configName, \
	HEAD_IDS(VT, DM, FT), \
	HEAD_INDEX(VT, DM, ALGO, FT), \
	VECTORS(VT, FT), \
	SSD_INDEX(VT, DM, ALGO, FT), \
	BUILD_SSD_BUILDER_CONFIG(VT, DM, ALGO), \
	SPTAG::VectorFileType::FT, \
	VECTOR_NUM, \
	VECTOR_DIM, \
	BS_internalResultNum, \
	BS_numberOfThreads, \
	BS_replicaCount, \
	BS_postingPageLimit, \
	true \
);} \

// DEFAULT
BDSSD(Float, L2, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Float, L2, KDT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Float, Cosine, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Float, Cosine, KDT, DEFAULT,
	60, 1, 4, 2)

BDSSD(Int8, L2, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int8, L2, KDT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, KDT, DEFAULT,
	60, 1, 4, 2)

BDSSD(UInt8, L2, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(UInt8, L2, KDT, DEFAULT,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, KDT, DEFAULT,
	60, 1, 4, 2)

BDSSD(Int16, L2, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int16, L2, KDT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, BKT, DEFAULT,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, KDT, DEFAULT,
	60, 1, 4, 2)


// XVEC
BDSSD(Float, L2, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Float, L2, KDT, XVEC,
	60, 1, 4, 2)
BDSSD(Float, Cosine, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Float, Cosine, KDT, XVEC,
	60, 1, 4, 2)

BDSSD(Int8, L2, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Int8, L2, KDT, XVEC,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, KDT, XVEC,
	60, 1, 4, 2)

BDSSD(UInt8, L2, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(UInt8, L2, KDT, XVEC,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, KDT, XVEC,
	60, 1, 4, 2)

BDSSD(Int16, L2, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Int16, L2, KDT, XVEC,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, BKT, XVEC,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, KDT, XVEC,
	60, 1, 4, 2)
#undef BDSSD


#define SCSSD(VT, DM, ALGO, FT, TFT, \
SS_internalResultNum, SS_resultNum) \
BOOST_AUTO_TEST_CASE(TestSearchSSDIndex##VT##DM##ALGO##FT##TFT) { \
string configName = SEARCH_SSD_CONFIG(VT, DM, ALGO); \
string queryFileName = QUERIES(VT, FT); \
string truthFileName = TRUTHSET(VT, DM, FT, TFT); \
string resultFileName = SEARCH_SSD_RESULT(VT, DM, ALGO, FT, TFT); \
TestSearchSSDIndex( \
	configName, \
	HEAD_IDS(VT, DM, FT), HEAD_INDEX(VT, DM, ALGO, FT), SSD_INDEX(VT, DM, ALGO, FT), SEARCH_SSD_BUILDER_CONFIG(VT, DM, ALGO), queryFileName, \
	SS_internalResultNum, SS_resultNum, 16, \
	resultFileName, \
	"",  \
	"", \
	truthFileName, \
	queryFileName, \
	"", \
	"10240", \
	0, \
	100000, \
	SPTAG::VectorFileType::FT, QUERY_NUM, VECTOR_DIM, \
	SPTAG::VectorFileType::FT, QUERY_NUM, VECTOR_DIM, \
	SPTAG::TruthFileType::TFT, QUERY_NUM \
);} \

SCSSD(Float, L2, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Float, L2, KDT, DEFAULT, TXT, 64, 64)
SCSSD(Float, Cosine, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Float, Cosine, KDT, DEFAULT, TXT, 64, 64)

SCSSD(Int8, L2, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Int8, L2, KDT, DEFAULT, TXT, 64, 64)
SCSSD(Int8, Cosine, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Int8, Cosine, KDT, DEFAULT, TXT, 64, 64)

SCSSD(UInt8, L2, BKT, DEFAULT, TXT, 64, 64)
SCSSD(UInt8, L2, KDT, DEFAULT, TXT, 64, 64)
SCSSD(UInt8, Cosine, BKT, DEFAULT, TXT, 64, 64)
SCSSD(UInt8, Cosine, KDT, DEFAULT, TXT, 64, 64)

SCSSD(Int16, L2, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Int16, L2, KDT, DEFAULT, TXT, 64, 64)
SCSSD(Int16, Cosine, BKT, DEFAULT, TXT, 64, 64)
SCSSD(Int16, Cosine, KDT, DEFAULT, TXT, 64, 64)


//Another
SCSSD(Float, L2, BKT, XVEC, XVEC, 64, 64)
SCSSD(Float, L2, KDT, XVEC, XVEC, 64, 64)
SCSSD(Float, Cosine, BKT, XVEC, XVEC, 64, 64)
SCSSD(Float, Cosine, KDT, XVEC, XVEC, 64, 64)

SCSSD(Int8, L2, BKT, XVEC, XVEC, 64, 64)
SCSSD(Int8, L2, KDT, XVEC, XVEC, 64, 64)
SCSSD(Int8, Cosine, BKT, XVEC, XVEC, 64, 64)
SCSSD(Int8, Cosine, KDT, XVEC, XVEC, 64, 64)

SCSSD(UInt8, L2, BKT, XVEC, XVEC, 64, 64)
SCSSD(UInt8, L2, KDT, XVEC, XVEC, 64, 64)
SCSSD(UInt8, Cosine, BKT, XVEC, XVEC, 64, 64)
SCSSD(UInt8, Cosine, KDT, XVEC, XVEC, 64, 64)

SCSSD(Int16, L2, BKT, XVEC, XVEC, 64, 64)
SCSSD(Int16, L2, KDT, XVEC, XVEC, 64, 64)
SCSSD(Int16, Cosine, BKT, XVEC, XVEC, 64, 64)
SCSSD(Int16, Cosine, KDT, XVEC, XVEC, 64, 64)
#undef SCSSD

BOOST_AUTO_TEST_SUITE_END()