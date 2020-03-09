#pragma once
#include "inc/SSDServing/Common/stdafx.h"
#include <cstring>
#include <limits>

using namespace std;

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {

    template <typename T>
    class CountVector
    {
        size_t m_bytes;
        T* m_data;
        T m_count;

    public:
        void Init(int size)
        {
            m_bytes = sizeof(T) * size;
            m_data = new T[size];
            m_count = 0;

            memset(m_data, 0, m_bytes);
        }

        CountVector() :m_data(nullptr) {}
        CountVector(int size) { Init(size); }
        ~CountVector() { if (m_data != nullptr) delete[] m_data; }

        inline void Clear()
        {
#undef max
            if (m_count == numeric_limits<T>::max())
            {
                memset(m_data, 0, m_bytes);
                m_count = 1;
            }
            else
            {
                m_count++;
            }
        }

        inline bool CheckAndSet(int idx)
        {
            if (m_data[idx] == m_count) return true;
            m_data[idx] = m_count;
            return false;
        }
    };

}
}
}