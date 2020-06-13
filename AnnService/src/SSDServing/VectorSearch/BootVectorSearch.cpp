#include "inc/SSDServing/VectorSearch/BootVectorSearch.h"
#include "inc/SSDServing/VectorSearch/BuildSsdIndex.h"
#include "inc/SSDServing/VectorSearch/SearchSsdIndex.h"
#include "inc/Helper/SimpleIniReader.h"
#include <iostream>

namespace SPTAG {
	namespace SSDServing {
		namespace VectorSearch {
			
			ErrorCode Bootstrap(Options& opts) {
                if (opts.m_buildSsdIndex)
                {
					std::cerr << "Start building SSD Index." << std::endl;
					if (false) {}
#define DefineVectorValueType(Name, Type) \
else if (COMMON_OPTS.m_valueType == VectorValueType::Name) { \
BuildSsdIndex<Type>(opts); \
} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
				}
				else {
					std::cerr << "Start searching SSD Index." << std::endl;
					if (false) {}
#define DefineVectorValueType(Name, Type) \
else if (COMMON_OPTS.m_valueType == VectorValueType::Name) { \
Search<Type>(opts); \
} \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType
				}
				
				return ErrorCode::Success;
			}
		}
	}
}