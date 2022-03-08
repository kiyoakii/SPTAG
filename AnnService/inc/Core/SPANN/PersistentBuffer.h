#include "inc/Helper/KeyValueIO.h"
#include <atomic>

namespace SPTAG {
    namespace SPANN {
        // concurrently safe with RocksDBIO
        class PersistentBuffer
        {
        public:
            PersistentBuffer(std::string& fileName, std::shared_ptr<Helper::KeyValueIO> db) : db(db), _size(0)
            {
                db->Initialize(fileName.c_str());
            }

            ~PersistentBuffer() {}

            inline int GetNewAssignmentID() { return _size++; }

            inline void GetAssignment(int assignmentID, std::string* assignment) { db->Get(assignmentID, assignment); }

            inline int PutAssignment(std::string& assignment)
            {
                int assignmentID = GetNewAssignmentID();
                db->Put(assignmentID, assignment);
                return assignmentID;
            }

            inline int GetCurrentAssignmentID() { return _size; }

        private:
            std::shared_ptr<Helper::KeyValueIO> db;
            std::atomic_int _size;
        };
    }
}