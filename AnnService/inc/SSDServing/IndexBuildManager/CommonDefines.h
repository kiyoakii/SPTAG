#pragma once

#include "Options.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"

namespace SPTAG {
	namespace SSDServing {
		extern BaseOptions COMMON_OPTS;
        using namespace ROCKSDB_NAMESPACE;
        extern std::string kDBPath;
        extern DB* db;
        extern Options dbOptions;
	}
}
