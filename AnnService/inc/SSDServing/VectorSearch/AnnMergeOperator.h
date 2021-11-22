#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "rocksdb/merge_operator.h"

namespace SPTAG {
	namespace SSDServing {

		class AnnMergeOperator : public ROCKSDB_NAMESPACE::AssociativeMergeOperator
        {
            public:
                virtual bool Merge(const ROCKSDB_NAMESPACE::Slice& key, const ROCKSDB_NAMESPACE::Slice* existing_value,
                     const ROCKSDB_NAMESPACE::Slice& value, std::string* new_value,
                     ROCKSDB_NAMESPACE::Logger* logger) const override {
                        std::string newPosting;
                        if(existing_value) 
                        {
                            newPosting += (*existing_value).ToString();
                            newPosting += value.ToString();
                        } else
                        {
                            newPosting += value.ToString();
                        }
                        *new_value = newPosting;
                        return true;
                    }
                virtual const char* Name() const override {
                    return "AnnMergeOperator";
                }
        };
    }
}