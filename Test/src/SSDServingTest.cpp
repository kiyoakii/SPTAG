#include "inc/Test.h"
#include <string>
#include <fstream>
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

#pragma omp parallel for
	for (int i = 0; i < querys.size(); ++i)
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
				exit(-1);
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
			exit(-1);
		}

		std::sort(neighbours.begin(), neighbours.end());

		for (size_t k = 0; k < K; k++)
		{
			truthset[i][k] = neighbours[k].key;
		}

	}
	
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

std::string CreateBaseConfig(SPTAG::VectorValueType p_valueType, SPTAG::DistCalcMethod p_distCalcMethod, SPTAG::DimensionType p_dim,
	std::string p_vectorPath, SPTAG::VectorFileType p_vectorType, SPTAG::SizeType p_vectorSize, std::string p_vectorDelimiter, 
	std::string p_queryPath, SPTAG::VectorFileType p_queryType, SPTAG::SizeType p_querySize, std::string p_queryDelimiter,
	std::string p_warmupPath, SPTAG::VectorFileType p_warmupType, SPTAG::SizeType p_warmupSize, std::string p_warmupDelimiter,
	std::string p_truthPath, SPTAG::TruthFileType p_truthType) {
	std::ostringstream config;
	config << "[Base]" << endl;
	config << "ValueType=" << SPTAG::Helper::Convert::ConvertToString(p_valueType) << endl;
	config << "DistCalcMethod=" << SPTAG::Helper::Convert::ConvertToString(p_distCalcMethod) << endl;
	config << "Dim=" << p_dim << endl;
	config << "VectorPath=" << p_vectorPath << endl;
	config << "VectorType=" << SPTAG::Helper::Convert::ConvertToString(p_vectorType) << endl;
	config << "VectorSize=" << p_vectorSize << endl;
	config << "VectorDelimiter=" << p_vectorDelimiter << endl;
	config << "QueryPath=" << p_queryPath << endl;
	config << "QueryType=" << SPTAG::Helper::Convert::ConvertToString(p_queryType) << endl;
	config << "QuerySize=" << p_querySize << endl;
	config << "QueryDelimiter=" << p_queryDelimiter << endl;
	config << "WarmupPath=" << p_warmupPath << endl;
	config << "WarmupType=" << SPTAG::Helper::Convert::ConvertToString(p_warmupType) << endl;
	config << "WarmupSize=" << p_warmupSize << endl;
	config << "WarmupDelimiter=" << p_warmupDelimiter << endl;
	config << "TruthPath=" << p_truthPath << endl;
	config << "TruthType=" << SPTAG::Helper::Convert::ConvertToString(p_truthType) << endl;
	config << endl;
	return config.str();
}

void TestHead(string configName, string OutputIDFile, string OutputVectorFile, string baseConfig) {

	ofstream config(configName);
	if (!config.is_open())
	{
		fprintf(stderr, "%s can't be opened. ", configName.c_str());
		BOOST_CHECK(false);
		return;
	}

	config << baseConfig;

	config << "[SelectHead]" << endl;
	config << "TreeNumber=" << "1" << endl;
	config << "BKTKmeansK=" << 3 << endl;
	config << "BKTLeafSize=" << 6 << endl;
	config << "SamplesNumber=" << 100 << endl;
	config << "NumberOfThreads=" << "2" << endl;
	config << "SaveBKT=" << "true" <<endl;

	config << "AnalyzeOnly=" << "false" << endl;
	config << "CalcStd=" << "true" << endl;
	config << "SelectDynamically=" << "true" << endl;
	config << "NoOutput=" << "false" << endl;

	config << "SelectThreshold=" << 12 << endl;
	config << "SplitFactor=" << 9 << endl;
	config << "SplitThreshold=" << 18 << endl;
	config << "Ratio=" << "0.2" << endl;
	config << "RecursiveCheckSmallCluster=" << "true" << endl;
	config << "PrintSizeCount=" << "true" << endl;

	config << "OutputIDFile=" << OutputIDFile << endl;
	config << "OutputVectorFile=" << OutputVectorFile << endl;

	config.close();

	char arg1[255], arg2[255];
	strncpy(arg1, "SSDServing", 255);
	strncpy(arg2, configName.c_str(), 255);
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
}

void TestBuildHead(
	string configName, 
	string p_headVectorFile, 
	string p_headIndexFile,
	SPTAG::IndexAlgoType p_indexAlgoType, 
	string p_builderFile, 
	string baseConfig) {

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
			bf.close();
		}
	}

	config << baseConfig;

	config << "[BuildHead]" << endl;
	config << "HeadVectorFile=" << p_headVectorFile << endl;
	config << "HeadIndex=" << p_headIndexFile << endl;
	config << "IndexAlgoType=" << SPTAG::Helper::Convert::ConvertToString(p_indexAlgoType) << endl;
	config << "BuilderConfigFile=" << p_builderFile << endl;
	config << "NumberOfThreads=" << 2 << endl;

	config.close();

	char arg1[255], arg2[255];
	strncpy(arg1, "SSDServing", 255);
	strncpy(arg2, configName.c_str(), 255);
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
}

void TestBuildSSDIndex(string configName,
	string p_vectorIDTranslate,
	string p_headIndexFolder,
	string p_headConfig,
	string p_ssdIndex,
	bool p_outputEmptyReplicaID,
	string baseConfig
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

	config << baseConfig;

	config << "[BuildSSDIndex]" << endl;
	config << "BuildSsdIndex=" << "true" << endl;
	config << "VectorIDTranslate=" << p_vectorIDTranslate << endl;
	config << "HeadIndexFolder=" << p_headIndexFolder << endl;
	config << "InternalResultNum=" << 60 << endl;
	config << "NumberOfThreads=" << 2 << endl;
	config << "HeadConfig=" << p_headConfig << endl;
	config << "SsdIndex=" << p_ssdIndex << endl;
	config << "ReplicaCount=" << 4 << endl;
	config << "PostingPageLimit=" << 2 << endl;
	config << "OutputEmptyReplicaID=" << p_outputEmptyReplicaID << endl;

	config.close();

	char arg1[255], arg2[255];
	strncpy(arg1, "SSDServing", 255);
	strncpy(arg2, configName.c_str(), 255);
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
}

void TestSearchSSDIndex(
	string configName,
	string p_vectorIDTranslate, 
	string p_headIndexFolder, 
	string p_headConfig,
	string p_searchResult,
	string p_extraFullGraphFile,
	string p_logFile,
	string baseConfig
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

	config << baseConfig;

	config << "[SearchSSDIndex]" << endl;
	config << "BuildSsdIndex=" << "false" << endl;
	config << "VectorIDTranslate=" << p_vectorIDTranslate << endl;
	config << "HeadIndexFolder=" << p_headIndexFolder << endl;
	config << "InternalResultNum=" << 64 << endl;
	config << "NumberOfThreads=" << 2 << endl;
	config << "HeadConfig=" << p_headConfig << endl;

	config << "SearchResult=" << p_searchResult << endl;
	config << "ExtraFullGraphFile=" << p_extraFullGraphFile << endl;
	config << "LogFile=" << p_logFile << endl;
	config << "QpsLimit=" << 0 << endl;
	config << "ResultNum=" << 64 << endl;
	config << "QueryCountLimit=" << 10000 << endl;

	config.close();

	char arg1[255], arg2[255];
	strncpy(arg1, "SSDServing", 255);
	strncpy(arg2, configName.c_str(), 255);
	char* params[2] = { arg1, arg2 };
	SPTAG::SSDServing::internalMain(2, params);
}

BOOST_AUTO_TEST_SUITE(SSDServingTest)

#define SSDTEST_DIRECTORY_NAME "sddtest"
#define VECTOR_NUM 1000
#define QUERY_NUM 10
#define VECTOR_DIM 100
//#if defined(_WIN32)
//#define SSDTEST_DIRECTORY SSDTEST_DIRECTORY_NAME "\\"
//#else
//#define SSDTEST_DIRECTORY SSDTEST_DIRECTORY_NAME "/"
//#endif
#define SSDTEST_DIRECTORY SSDTEST_DIRECTORY_NAME "/"

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

#define WTEV(VT, DM, FT) \
BOOST_AUTO_TEST_CASE(TestHead##VT##DM##FT) \
{ \
	string configName = SELECT_HEAD_CONFIG(VT, DM, FT); \
	string OutputIDFile = HEAD_IDS(VT, DM, FT); \
	string OutputVectorFile = HEAD_VECTORS(VT, DM, FT); \
	string base_config = CreateBaseConfig(SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, VECTOR_DIM, \
		VECTORS(VT, FT), SPTAG::VectorFileType::FT, VECTOR_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		"", SPTAG::TruthFileType::Undefined \
		); \
\
	TestHead(configName, OutputIDFile, OutputVectorFile, base_config);  \
} \

WTEV(Float, L2, DEFAULT)
WTEV(Float, Cosine, DEFAULT)
WTEV(Int16, L2, DEFAULT)
WTEV(Int16, Cosine, DEFAULT)
WTEV(UInt8, L2, DEFAULT)
WTEV(UInt8, Cosine, DEFAULT)
WTEV(Int8, L2, DEFAULT)
WTEV(Int8, Cosine, DEFAULT)

WTEV(Float, L2, XVEC)
WTEV(Float, Cosine, XVEC)
WTEV(Int16, L2, XVEC)
WTEV(Int16, Cosine, XVEC)
WTEV(UInt8, L2, XVEC)
WTEV(UInt8, Cosine, XVEC)
WTEV(Int8, L2, XVEC)
WTEV(Int8, Cosine, XVEC)

#undef WTEV

#define BDHD(VT, DM, ALGO, FT) \
BOOST_AUTO_TEST_CASE(TestBuildHead##VT##DM##ALGO##FT) { \
	string configName = BUILD_HEAD_CONFIG(VT, DM, ALGO); \
	string builderFile = BUILD_HEAD_BUILDER_CONFIG(VT, DM, ALGO); \
	string base_config = CreateBaseConfig(SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, VECTOR_DIM, \
		VECTORS(VT, FT), SPTAG::VectorFileType::FT, VECTOR_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		"", SPTAG::TruthFileType::Undefined \
		); \
TestBuildHead( \
	configName, \
	HEAD_VECTORS(VT, DM, FT), \
	HEAD_INDEX(VT, DM, ALGO, FT), \
	SPTAG::IndexAlgoType::ALGO, \
	builderFile, \
	base_config \
);} \

BDHD(Float, L2, BKT, DEFAULT)
BDHD(Float, L2, KDT, DEFAULT)
BDHD(Float, Cosine, BKT, DEFAULT)
BDHD(Float, Cosine, KDT, DEFAULT)

BDHD(Int8, L2, BKT, DEFAULT)
BDHD(Int8, L2, KDT, DEFAULT)
BDHD(Int8, Cosine, BKT, DEFAULT)
BDHD(Int8, Cosine, KDT, DEFAULT)

BDHD(UInt8, L2, BKT, DEFAULT)
BDHD(UInt8, L2, KDT, DEFAULT)
BDHD(UInt8, Cosine, BKT, DEFAULT)
BDHD(UInt8, Cosine, KDT, DEFAULT)

BDHD(Int16, L2, BKT, DEFAULT)
BDHD(Int16, L2, KDT, DEFAULT)
BDHD(Int16, Cosine, BKT, DEFAULT)
BDHD(Int16, Cosine, KDT, DEFAULT)

//XVEC
BDHD(Float, L2, BKT, XVEC)
BDHD(Float, L2, KDT, XVEC)
BDHD(Float, Cosine, BKT, XVEC)
BDHD(Float, Cosine, KDT, XVEC)

BDHD(Int8, L2, BKT, XVEC)
BDHD(Int8, L2, KDT, XVEC)
BDHD(Int8, Cosine, BKT, XVEC)
BDHD(Int8, Cosine, KDT, XVEC)

BDHD(UInt8, L2, BKT, XVEC)
BDHD(UInt8, L2, KDT, XVEC)
BDHD(UInt8, Cosine, BKT, XVEC)
BDHD(UInt8, Cosine, KDT, XVEC)

BDHD(Int16, L2, BKT, XVEC)
BDHD(Int16, L2, KDT, XVEC)
BDHD(Int16, Cosine, BKT, XVEC)
BDHD(Int16, Cosine, KDT, XVEC)

#undef BDHD

#define BDSSD(VT, DM, ALGO, FT) \
BOOST_AUTO_TEST_CASE(TestBuildSSDIndex##VT##DM##ALGO##FT) { \
	string configName = BUILD_SSD_CONFIG(VT, DM, ALGO); \
	string base_config = CreateBaseConfig(SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, VECTOR_DIM, \
		VECTORS(VT, FT), SPTAG::VectorFileType::FT, VECTOR_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		"", SPTAG::TruthFileType::Undefined \
		); \
TestBuildSSDIndex(\
	configName, \
	HEAD_IDS(VT, DM, FT), \
	HEAD_INDEX(VT, DM, ALGO, FT), \
	BUILD_SSD_BUILDER_CONFIG(VT, DM, ALGO), \
	SSD_INDEX(VT, DM, ALGO, FT), \
	true, \
	base_config \
);} \

// DEFAULT
BDSSD(Float, L2, BKT, DEFAULT)
BDSSD(Float, L2, KDT, DEFAULT)
BDSSD(Float, Cosine, BKT, DEFAULT)
BDSSD(Float, Cosine, KDT, DEFAULT)

BDSSD(Int8, L2, BKT, DEFAULT)
BDSSD(Int8, L2, KDT, DEFAULT)
BDSSD(Int8, Cosine, BKT, DEFAULT)
BDSSD(Int8, Cosine, KDT, DEFAULT)

BDSSD(UInt8, L2, BKT, DEFAULT)
BDSSD(UInt8, L2, KDT, DEFAULT)
BDSSD(UInt8, Cosine, BKT, DEFAULT)
BDSSD(UInt8, Cosine, KDT, DEFAULT)

BDSSD(Int16, L2, BKT, DEFAULT)
BDSSD(Int16, L2, KDT, DEFAULT)
BDSSD(Int16, Cosine, BKT, DEFAULT)
BDSSD(Int16, Cosine, KDT, DEFAULT)

// XVEC
BDSSD(Float, L2, BKT, XVEC)
BDSSD(Float, L2, KDT, XVEC)
BDSSD(Float, Cosine, BKT, XVEC)
BDSSD(Float, Cosine, KDT, XVEC)

BDSSD(Int8, L2, BKT, XVEC)
BDSSD(Int8, L2, KDT, XVEC)
BDSSD(Int8, Cosine, BKT, XVEC)
BDSSD(Int8, Cosine, KDT, XVEC)

BDSSD(UInt8, L2, BKT, XVEC)
BDSSD(UInt8, L2, KDT, XVEC)
BDSSD(UInt8, Cosine, BKT, XVEC)
BDSSD(UInt8, Cosine, KDT, XVEC)

BDSSD(Int16, L2, BKT, XVEC)
BDSSD(Int16, L2, KDT, XVEC)
BDSSD(Int16, Cosine, BKT, XVEC)
BDSSD(Int16, Cosine, KDT, XVEC)
#undef BDSSD


#define SCSSD(VT, DM, ALGO, FT, TFT) \
BOOST_AUTO_TEST_CASE(TestSearchSSDIndex##VT##DM##ALGO##FT##TFT) { \
	string configName = SEARCH_SSD_CONFIG(VT, DM, ALGO); \
	string base_config = CreateBaseConfig(SPTAG::VectorValueType::VT, SPTAG::DistCalcMethod::DM, VECTOR_DIM, \
		VECTORS(VT, FT), SPTAG::VectorFileType::FT, VECTOR_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		QUERIES(VT, FT), SPTAG::VectorFileType::FT, QUERY_NUM, "", \
		TRUTHSET(VT, DM, FT, TFT), SPTAG::TruthFileType::TFT \
		); \
\
TestSearchSSDIndex( \
	configName, \
	HEAD_IDS(VT, DM, FT), \
	HEAD_INDEX(VT, DM, ALGO, FT), \
	SEARCH_SSD_BUILDER_CONFIG(VT, DM, ALGO), \
	SEARCH_SSD_RESULT(VT, DM, ALGO, FT, TFT), \
	SSD_INDEX(VT, DM, ALGO, FT), \
	"", \
	base_config \
);} \

SCSSD(Float, L2, BKT, DEFAULT, TXT)
SCSSD(Float, L2, KDT, DEFAULT, TXT)
SCSSD(Float, Cosine, BKT, DEFAULT, TXT)
SCSSD(Float, Cosine, KDT, DEFAULT, TXT)

SCSSD(Int8, L2, BKT, DEFAULT, TXT)
SCSSD(Int8, L2, KDT, DEFAULT, TXT)
SCSSD(Int8, Cosine, BKT, DEFAULT, TXT)
SCSSD(Int8, Cosine, KDT, DEFAULT, TXT)

SCSSD(UInt8, L2, BKT, DEFAULT, TXT)
SCSSD(UInt8, L2, KDT, DEFAULT, TXT)
SCSSD(UInt8, Cosine, BKT, DEFAULT, TXT)
SCSSD(UInt8, Cosine, KDT, DEFAULT, TXT)

SCSSD(Int16, L2, BKT, DEFAULT, TXT)
SCSSD(Int16, L2, KDT, DEFAULT, TXT)
SCSSD(Int16, Cosine, BKT, DEFAULT, TXT)
SCSSD(Int16, Cosine, KDT, DEFAULT, TXT)


//Another
SCSSD(Float, L2, BKT, XVEC, XVEC)
SCSSD(Float, L2, KDT, XVEC, XVEC)
SCSSD(Float, Cosine, BKT, XVEC, XVEC)
SCSSD(Float, Cosine, KDT, XVEC, XVEC)

SCSSD(Int8, L2, BKT, XVEC, XVEC)
SCSSD(Int8, L2, KDT, XVEC, XVEC)
SCSSD(Int8, Cosine, BKT, XVEC, XVEC)
SCSSD(Int8, Cosine, KDT, XVEC, XVEC)

SCSSD(UInt8, L2, BKT, XVEC, XVEC)
SCSSD(UInt8, L2, KDT, XVEC, XVEC)
SCSSD(UInt8, Cosine, BKT, XVEC, XVEC)
SCSSD(UInt8, Cosine, KDT, XVEC, XVEC)

SCSSD(Int16, L2, BKT, XVEC, XVEC)
SCSSD(Int16, L2, KDT, XVEC, XVEC)
SCSSD(Int16, Cosine, BKT, XVEC, XVEC)
SCSSD(Int16, Cosine, KDT, XVEC, XVEC)
#undef SCSSD

BOOST_AUTO_TEST_SUITE_END()