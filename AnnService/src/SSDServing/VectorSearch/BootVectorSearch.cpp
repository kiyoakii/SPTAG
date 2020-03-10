#include "inc/SSDServing/Common/stdafx.h"
#include "inc/SSDServing/VectorSearch/BootVectorSearch.h"
#include "inc/SSDServing/VectorSearch/BuildSsdIndex.h"
#include "inc/SSDServing/VectorSearch/SearchSsdIndex.h"
#include "inc/Helper/SimpleIniReader.h"

namespace SPTAG {
	namespace SSDServing {
		namespace VectorSearch {

			void GetHeadIndex(Options& p_opts, shared_ptr<VectorIndex>& p_index) {
				if (VectorIndex::LoadIndex(p_opts.m_headIndexFolder, p_index) != ErrorCode::Success) {
					std::cerr << "ERROR: Cannot Load index files!" << std::endl;
					exit(1);
				}
				p_index->SetParameter("NumberOfThreads", std::to_string(p_opts.m_iNumberOfThreads));
				Helper::IniReader iniReader;
				if (!p_opts.m_headConfig.empty())
				{
					if (iniReader.LoadIniFile(p_opts.m_headConfig) != ErrorCode::Success) {
						std::cerr << "ERROR of loading head index config: " << p_opts.m_headConfig << std::endl;
						exit(1);
					}

					for (const auto& iter : iniReader.GetParameters("Index"))
					{
						p_index->SetParameter(iter.first.c_str(), iter.second.c_str());
					}
				}
			}

			ErrorCode Bootstrap(Options& opts) {
				shared_ptr<VectorIndex> headIndex;
				GetHeadIndex(opts, headIndex);

                if (opts.m_buildSsdIndex)
                {
					cerr << "Start building SSD Index." << endl;
					if (false) {}
#define DefineVectorValueType(Name, Type) \
else if (headIndex->GetVectorValueType() == VectorValueType::Name) { \
BuildSsdIndex<Type>(opts, headIndex); \
} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
				}
				else {
					cerr << "Start searching SSD Index." << endl;
					if (false) {}
#define DefineVectorValueType(Name, Type) \
else if (headIndex->GetVectorValueType() == VectorValueType::Name) { \
Search<Type>(opts, headIndex); \
} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
				}
				
				return ErrorCode::Success;
			}
		}
	}
}