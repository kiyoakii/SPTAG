#include "inc/Helper/DiskIO.h"
#include <atomic>

namespace SPTAG {
    namespace SPANN {
        class PersistentBuffer
        {
        public:
            PersistentBuffer(std::string& fileName)
            {
                db->Initialize(fileName.c_str());
                m_updateID.store(0);
            }

            ~PersistentBuffer() {}

            int GetNewAssignmentID()
            {
                return m_updateID.fetch_add(1);
            }

            void GetAssignment(int assignmentID, std::string* assignment)
            {
                db->Get(assignmentID, assignment);
            }

            int PutAssignment(std::string& assignment)
            {
                int assignmentID = GetNewAssignmentID();
                db->Put(assignmentID, assignment);
                return assignmentID;
            }

            int GetCurrentAssignmentID()
            {
                return m_updateID.load();
            }

        private:
            std::shared_ptr<Helper::KeyValueIO> db;
            //update assignment id
            std::atomic_int m_updateID;
        };
    }
}