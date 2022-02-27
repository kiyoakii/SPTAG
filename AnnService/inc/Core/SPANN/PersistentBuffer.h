#include "assert.h"
#include "rocksdb/db.h"
#include "rocksdb/slice.h"
#include "rocksdb/options.h"
#include "inc/Helper/StringConvert.h"

namespace SPTAG {
    namespace PANNS {
        namespace FRESH {
            class PersistentBuffer
            {
                public:
                    PersistentBuffer(std::string& fileName)
                    {
                        db.Initialize(fileName.c_str());
                        m_updateID.store(0);
                    }

                    ~PersistentBuffer() {}

                    int GetNewAssignmentID()
                    {
                        return m_updateID.fetch_add(1);
                    }

                    void GetAssignment(int assignmentID, std::string* assignment)
                    {
                        db.Get(assignmentID, assignment);
                    }

                    int PutAssignment(std::string& assignment)
                    {
                        int assignmentID = GetNewAssignmentID();
                        db.Put(assignmentID, assignment);
                        return assignmentID;
                    }

                    int GetCurrentAssignmentID()
                    {
                        return m_updateID.load();
                    }

                private:
                    Helper::RocksDBIO db;
                    //update assignment id
                    std::atomic_int m_updateID;
            };
        }
    }
}