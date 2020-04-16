#pragma once

#include <cstring>
#include <limits>

using namespace std;

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {

            class HashBasedDeduper
            {
            public:
                HashBasedDeduper()
                {
                    Init();
                }

                ~HashBasedDeduper() {}

                void Init(int stub = 0)
                {
                    m_secondHash = true;
                    Clear();
                }

                void Clear()
                {
                    if (!m_secondHash)
                    {
                        // Clear first block.
                        memset(&m_hashTable[0], 0, sizeof(int) * (m_poolSize + 1));
                    }
                    else
                    {
                        // Clear all blocks.
                        memset(&m_hashTable[0], 0, 2 * sizeof(int) * (m_poolSize + 1));
                        m_secondHash = false;
                    }
                }

                inline bool CheckAndSet(int idx)
                {
                    // Inner Index is begin from 1
                    return _CheckAndSet(&m_hashTable[0], idx + 1) == 0;
                }

            protected:
                // Max loop number in one hash block.
                static const int m_maxLoop = 8;

                // Max pool size.
                static const int m_poolSize = 8191;

                // Could we use the second hash block.
                bool m_secondHash;

                // Record 2 hash tables.
                // [0~m_poolSize + 1) is the first block.
                // [m_poolSize + 1, 2*(m_poolSize + 1)) is the second block;
                int m_hashTable[(m_poolSize + 1) * 2];


                inline unsigned hash_func2(int idx, int loop)
                {
                    return ((unsigned)idx + loop) & m_poolSize;
                }


                inline unsigned hash_func(unsigned idx)
                {
                    return ((unsigned)(idx * 99991) + _rotl(idx, 2) + 101) & m_poolSize;
                }

                inline int _CheckAndSet(int* hashTable, int idx)
                {
                    unsigned index, loop;

                    // Get first hash position.
                    index = hash_func(idx);
                    for (loop = 0; loop < m_maxLoop; ++loop)
                    {
                        if (!hashTable[index])
                        {
                            // index first match and record it.
                            hashTable[index] = idx;
                            return 1;
                        }
                        if (hashTable[index] == idx)
                        {
                            // Hit this item in hash table.
                            return 0;
                        }
                        // Get next hash position.
                        index = hash_func2(index, loop);
                    }

                    if (hashTable == &m_hashTable[0])
                    {
                        // Use second hash block.
                        m_secondHash = true;
                        return _CheckAndSet(&m_hashTable[m_poolSize + 1], idx);
                    }

                    // Do not include this item.
                    return -1;
                }
            };

}
}
}