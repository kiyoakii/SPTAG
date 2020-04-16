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
void GenerateVectors(string fileName, SPTAG::SizeType rows, SPTAG::DimensionType dims) {
	if (boost::filesystem::exists(fileName))
	{
		fprintf(stdout, "%s was generated. Skip generation.", fileName.c_str());
		return;
	}
	uniform_real_distribution<float> ud(0, 126);
	ofstream of(fileName, ofstream::binary);
	if (!of.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", fileName.c_str());
		BOOST_CHECK(false);
		return;
	}
	of.write(reinterpret_cast<char*>(&rows), sizeof(rows));
	of.write(reinterpret_cast<char*>(&dims), sizeof(dims));
	
	mt19937 mt(543);
	for (size_t i = 0; i < rows; i++)
	{
		vector<T> tmp(dims, 0);
		for (size_t j = 0; j < dims; j++)
		{
			float smt = ud(mt);
			tmp[j] = static_cast<T>(smt);
		}

		SPTAG::COMMON::Utils::Normalize(tmp.data(), dims, SPTAG::COMMON::Utils::GetBase<T>());
		of.write(reinterpret_cast<char*>(tmp.data()), tmp.size() * sizeof(T));
	}

	of.close();
}

void GenVec(string vectorsName, SPTAG::VectorValueType vecType, SPTAG::SizeType rows = 1000, SPTAG::DimensionType dims = 100) {
	switch (vecType)
	{
#define DefineVectorValueType(Name, Type) \
case SPTAG::VectorValueType::Name: \
GenerateVectors<Type>(vectorsName, rows, dims); \
break; \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
	default:
		break;
	}
}

template<typename T>
void FillVectors(const string binName, vector<vector<T>>& bins) {
	ifstream f(binName, ifstream::binary);
	if (!f.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", binName.c_str());
		BOOST_CHECK(false);
		return;
	}
	SPTAG::SizeType rowNum;
	SPTAG::DimensionType dimNum;
	f.read(reinterpret_cast<char*>(&rowNum), sizeof(SPTAG::SizeType));
	f.read(reinterpret_cast<char*>(&dimNum), sizeof(SPTAG::DimensionType));
	for (size_t i = 0; i < rowNum; i++)
	{
		vector<T> cur(dimNum, 0);
		f.read(reinterpret_cast<char*>(cur.data()), dimNum * sizeof(T));
		bins.push_back(cur);
	}
	f.close();
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

template<typename T>
void GenerateTruth(const string queryFile, const string vectorFile, const string truthFile,
	const SPTAG::DistCalcMethod distMethod, const int K) {
	if (boost::filesystem::exists(truthFile))
	{
		fprintf(stdout, "truthFile: %s was generated. Skip generation.", truthFile.c_str());
		return;
	}
	vector<vector<T>> querys;
	vector<vector<T>> vectors;
	FillVectors<T>(queryFile, querys);
	FillVectors<T>(vectorFile, vectors);
	vector<vector<SPTAG::SizeType>> truthset(querys.size(), vector<SPTAG::SizeType>(K, 0));
	ofstream of(truthFile);
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

	for (size_t i = 0; i < querys.size(); i++)
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

	of.close();
	
}

void GenerateTruth(const string queryFile, const string vectorFile, const string truthFile,
	const SPTAG::DistCalcMethod distMethod, const int K, SPTAG::VectorValueType vvt) {
#define DefineVectorValueType(Name, Type) \
	if (vvt == SPTAG::VectorValueType::Name) { \
		GenerateTruth<Type>(queryFile, vectorFile, truthFile, distMethod, K); \
	} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
}

void TestHead(string vectorsName, string configName, string OutputIDFile, string OutputVectorFile, 
	SPTAG::VectorValueType vecType, SPTAG::DistCalcMethod distMethod,
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

	config << "TreeNumber=" << "1" << endl;
	config << "BKTKmeansK=" << p_BKTKmeansK << endl;
	config << "BKTLeafSize=" << p_BKTLeafSize << endl;
	config << "SamplesNumber=" << p_SamplesNumber << endl;
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

	config << "SsdIndex=" << p_ssdIndex << endl;
	config << "HeadConfig=" << p_headConfig << endl;
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
	int p_queryCountLimit
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
#define SSDTEST_DIRECTORY SSDTEST_DIRECTORY_NAME "\\"
#define VECTORS(VT) SSDTEST_DIRECTORY "vectors_"#VT".bin"
#define QUERIES(VT) SSDTEST_DIRECTORY "vectors_"#VT".query"
#define TRUTHSET(VT, DM) SSDTEST_DIRECTORY "vectors_"#VT"_"#DM".truth"
#define HEAD_IDS(VT, DM) SSDTEST_DIRECTORY "head_ids_"#VT"_"#DM".bin"
#define HEAD_VECTORS(VT, DM) SSDTEST_DIRECTORY "head_vectors_"#VT"_"#DM".bin"
#define HEAD_INDEX(VT, DM, ALGO) SSDTEST_DIRECTORY "head_"#VT"_"#DM"_"#ALGO".head_index"
#define SSD_INDEX(VT, DM, ALGO) SSDTEST_DIRECTORY "ssd_"#VT"_"#DM"_"#ALGO".ssd_index"

#define SELECT_HEAD_CONFIG(VT, DM) SSDTEST_DIRECTORY "test_head_"#VT"_"#DM".ini"
#define BUILD_HEAD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_head_"#VT"_"#DM"_"#ALGO".ini"
#define BUILD_HEAD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_head_"#VT"_"#DM"_"#ALGO".builder.ini"
#define BUILD_SSD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_ssd"#VT"_"#DM"_"#ALGO".ini"
#define BUILD_SSD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_build_ssd"#VT"_"#DM"_"#ALGO".builder.ini"
#define SEARCH_SSD_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO".ini"
#define SEARCH_SSD_BUILDER_CONFIG(VT, DM, ALGO) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO".builder.ini"
#define SEARCH_SSD_RESULT(VT, DM, ALGO) SSDTEST_DIRECTORY "test_search_ssd_"#VT"_"#DM"_"#ALGO".result"

#define GVQ(VT) \
BOOST_AUTO_TEST_CASE(GenerateVectorsQueries##VT) { \
boost::filesystem::create_directory(SSDTEST_DIRECTORY_NAME); \
GenVec(VECTORS(VT), SPTAG::VectorValueType::VT, 1000, 100); \
GenVec(QUERIES(VT), SPTAG::VectorValueType::VT, 10, 100); \
} \

GVQ(Float)
GVQ(Int16)
GVQ(UInt8)
GVQ(Int8)
#undef GVQ

#define GTR(VT, DM) \
BOOST_AUTO_TEST_CASE(GenerateTruth##VT##DM) { \
GenerateTruth(QUERIES(VT), VECTORS(VT), TRUTHSET(VT, DM), \
	SPTAG::DistCalcMethod::DM, 128, SPTAG::VectorValueType::VT); \
} \

GTR(Float, L2)
GTR(Float, Cosine)
GTR(Int16, L2)
GTR(Int16, Cosine)
GTR(UInt8, L2)
GTR(UInt8, Cosine)
GTR(Int8, L2)
GTR(Int8, Cosine)
#undef GTR

#define WTEV(VT, DM, \
	S_BKTKmeansK, S_BKTLeafSize, S_SamplesNumber, \
S_SelectThreshold, S_SplitFactor, S_SplitThreshold) \
BOOST_AUTO_TEST_CASE(TestHead##VT##DM) \
{ \
	string vectorsName = VECTORS(VT); \
	string configName = SELECT_HEAD_CONFIG(VT, DM); \
	string OutputIDFile = HEAD_IDS(VT, DM); \
	string OutputVectorFile = HEAD_VECTORS(VT, DM); \
\
	TestHead(vectorsName, configName, OutputIDFile, OutputVectorFile, \
		SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, \
		S_BKTKmeansK, S_BKTLeafSize, S_SamplesNumber, \
		S_SelectThreshold, S_SplitFactor, S_SplitThreshold \
	); \
} \

WTEV(Float, L2, 
	3, 6, 100, 
	12, 9, 18)
WTEV(Float, Cosine,
	3, 6, 100, 
	12, 9, 18)
WTEV(Int16, L2,
	3, 6, 100,
	12, 9, 18)
WTEV(Int16, Cosine,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, L2,
	3, 6, 100,
	12, 9, 18)
WTEV(UInt8, Cosine,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, L2,
	3, 6, 100,
	12, 9, 18)
WTEV(Int8, Cosine,
	3, 6, 100,
	12, 9, 18)

#undef WTEV

#define BDHD(VT, DM, ALGO, \
H_BKTKmeansK, H_BKTLeafSize, H_Samples) \
BOOST_AUTO_TEST_CASE(TestBuildHead##VT##DM##ALGO) { \
string configName = BUILD_HEAD_CONFIG(VT, DM, ALGO); \
string builderFile = BUILD_HEAD_BUILDER_CONFIG(VT, DM, ALGO); \
TestBuildHead( \
	configName, \
	HEAD_VECTORS(VT, DM), \
	HEAD_INDEX(VT, DM, ALGO), \
	SPTAG::IndexAlgoType::ALGO, \
	SPTAG::DistCalcMethod::DM, \
	builderFile, \
	SPTAG::VectorValueType::VT, \
	H_BKTKmeansK, H_BKTLeafSize, H_Samples, \
	3 \
);} \

BDHD(Float, L2, BKT, 
	2, 4, 20)
BDHD(Float, L2, KDT,
	2, 4, 20)
BDHD(Float, Cosine, BKT,
	2, 4, 20)
BDHD(Float, Cosine, KDT,
	2, 4, 20)

BDHD(Int8, L2, BKT,
	2, 4, 20)
BDHD(Int8, L2, KDT,
	2, 4, 20)
BDHD(Int8, Cosine, BKT,
	2, 4, 20)
BDHD(Int8, Cosine, KDT,
	2, 4, 20)

BDHD(UInt8, L2, BKT,
	2, 4, 20)
BDHD(UInt8, L2, KDT,
	2, 4, 20)
BDHD(UInt8, Cosine, BKT,
	2, 4, 20)
BDHD(UInt8, Cosine, KDT,
	2, 4, 20)

BDHD(Int16, L2, BKT,
	2, 4, 20)
BDHD(Int16, L2, KDT,
	2, 4, 20)
BDHD(Int16, Cosine, BKT,
	2, 4, 20)
BDHD(Int16, Cosine, KDT,
	2, 4, 20)

#undef BDHD

#define BDSSD(VT, DM, ALGO, \
	BS_internalResultNum, BS_numberOfThreads, BS_replicaCount, BS_postingPageLimit) \
BOOST_AUTO_TEST_CASE(TestBuildSSDIndex##VT##DM##ALGO) { \
string configName = BUILD_SSD_CONFIG(VT, DM, ALGO); \
TestBuildSSDIndex(\
	configName, \
	HEAD_IDS(VT, DM), \
	HEAD_INDEX(VT, DM, ALGO), \
	VECTORS(VT), \
	SSD_INDEX(VT, DM, ALGO), \
	BUILD_SSD_BUILDER_CONFIG(VT, DM, ALGO), \
	BS_internalResultNum, \
	BS_numberOfThreads, \
	BS_replicaCount, \
	BS_postingPageLimit, \
	true \
);} \

BDSSD(Float, L2, BKT, 
	60, 1, 4, 2)
BDSSD(Float, L2, KDT,
	60, 1, 4, 2)
BDSSD(Float, Cosine, BKT,
	60, 1, 4, 2)
BDSSD(Float, Cosine, KDT,
	60, 1, 4, 2)

BDSSD(Int8, L2, BKT,
	60, 1, 4, 2)
BDSSD(Int8, L2, KDT,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, BKT,
	60, 1, 4, 2)
BDSSD(Int8, Cosine, KDT,
	60, 1, 4, 2)

BDSSD(UInt8, L2, BKT,
	60, 1, 4, 2)
BDSSD(UInt8, L2, KDT,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, BKT,
	60, 1, 4, 2)
BDSSD(UInt8, Cosine, KDT,
	60, 1, 4, 2)

BDSSD(Int16, L2, BKT,
	60, 1, 4, 2)
BDSSD(Int16, L2, KDT,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, BKT,
	60, 1, 4, 2)
BDSSD(Int16, Cosine, KDT,
	60, 1, 4, 2)
#undef BDSSD


#define SCSSD(VT, DM, ALGO, SS_internalResultNum, SS_resultNum) \
BOOST_AUTO_TEST_CASE(TestSearchSSDIndex##VT##DM##ALGO) { \
string configName = SEARCH_SSD_CONFIG(VT, DM, ALGO); \
string queryFileName = QUERIES(VT); \
string truthFileName = TRUTHSET(VT, DM); \
string resultFileName = SEARCH_SSD_RESULT(VT, DM, ALGO); \
TestSearchSSDIndex( \
	configName, \
	HEAD_IDS(VT, DM), HEAD_INDEX(VT, DM, ALGO), SSD_INDEX(VT, DM, ALGO), SEARCH_SSD_BUILDER_CONFIG(VT, DM, ALGO), queryFileName, \
	SS_internalResultNum, SS_resultNum, 16, \
	resultFileName, \
	"",  \
	"", \
	truthFileName, \
	"", \
	"", \
	"10240", \
	0, \
	100000 \
);} \

SCSSD(Float, L2, BKT, 64, 64)
SCSSD(Float, L2, KDT, 64, 64)
SCSSD(Float, Cosine, BKT, 64, 64)
SCSSD(Float, Cosine, KDT, 64, 64)

SCSSD(Int8, L2, BKT, 64, 64)
SCSSD(Int8, L2, KDT, 64, 64)
SCSSD(Int8, Cosine, BKT, 64, 64)
SCSSD(Int8, Cosine, KDT, 64, 64)

SCSSD(UInt8, L2, BKT, 64, 64)
SCSSD(UInt8, L2, KDT, 64, 64)
SCSSD(UInt8, Cosine, BKT, 64, 64)
SCSSD(UInt8, Cosine, KDT, 64, 64)

SCSSD(Int16, L2, BKT, 64, 64)
SCSSD(Int16, L2, KDT, 64, 64)
SCSSD(Int16, Cosine, BKT, 64, 64)
SCSSD(Int16, Cosine, KDT, 64, 64)
#undef SCSSD

BOOST_AUTO_TEST_SUITE_END()