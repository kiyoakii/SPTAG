#pragma once

#include "Options.h"
#include "inc/rocksdb/db.h"
#include "inc/rocksdb/slice.h"
#include "inc/rocksdb/options.h"

namespace SPTAG {
	namespace SSDServing {
		extern BaseOptions COMMON_OPTS;
        using namespace ROCKSDB_NAMESPACE;
        std::string kDBPath = "/tmp/rocksdb_simple_example";
        DB* db;
        Options dbOptions;
	}
}