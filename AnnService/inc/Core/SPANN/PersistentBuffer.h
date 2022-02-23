#include "assert.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "inc/Helper/StringConvert.h"

namespace SPTAG {
    namespace SPANN {
        class PersistentBuffer
        {
            public:
                PersistentBuffer(std::string& fileName)
                {
                    pDBpath = fileName;
                    pdbOptions.IncreaseParallelism();
                    pdbOptions.OptimizeLevelStyleCompaction();
                    pdbOptions.create_if_missing = true;
                    ROCKSDB_NAMESPACE::Status ps = ROCKSDB_NAMESPACE::DB::Open(pdbOptions, pDBpath, &pDB);
                    assert(ps.op());
                    m_updateID.store(0);
                }

                int GetNewAssignmentID()
                {
                    return m_updateID.fetch_add(1);
                }

                void GetAssignment(int assignmentID, std::string* assignment)
                {
                    pDB->Get(ReadOptions(), Helper::Convert::ConvertToString<int>(assignmentID), assignment);
                }

                int PutAssignment(std::string& assignment)
                {
                    int assignmentID = GetNewAssignmentID();
                    pDB->Put(WriteOptions(), Helper::Convert::ConvertToString<int>(assignmentID), assignment);
                    return assignmentID;
                }

                int GetCurrentAssignmentID()
                {
                    return m_updateID.load();
                }

                void StopPDB()
                {
                    pDB->Close();
                    ROCKSDB_NAMESPACE::DestroyDB(pDBpath, pdbOptions);
                    delete pDB;
                }

            private:
                ROCKSDB_NAMESPACE::DB* pDB;
                ROCKSDB_NAMESPACE::Options pdbOptions;
                std::string pDBpath;
                //update assignment id
                std::atomic_int m_updateID;
        };
    }
}