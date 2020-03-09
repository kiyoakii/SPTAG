#pragma once
#include "inc/SSDServing/Common/stdafx.h"
#include "inc/SSDServing/VectorSearch/SearchStats.h"
#include "inc/SSDServing/VectorSearch/VectorSearchUtils.h"
#include "inc/SSDServing/VectorSearch/Options.h"
#include "inc/Core/Common/QueryResultSet.h"

#include <memory>
#include <vector>
#include <functional>
#include <windows.h>

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch{
            template<typename T>
            class PageBuffer
            {
            public:
                PageBuffer()
                    : m_pageBufferSize(0)
                {
                }

                void ReservePageBuffer(std::size_t p_size)
                {
                    if (m_pageBufferSize < p_size)
                    {
                        m_pageBufferSize = p_size;
                        m_pageBuffer.reset(new T[m_pageBufferSize], std::default_delete<T[]>());
                    }
                }

                T* GetBuffer()
                {
                    return m_pageBuffer.get();
                }


            private:
                std::shared_ptr<T> m_pageBuffer;

                std::size_t m_pageBufferSize;
            };


            struct ExtraWorkSpace
            {
                ExtraWorkSpace()
                    : m_fileHandle(NULL), m_fileHandle2(NULL)
                {
                }

                ~ExtraWorkSpace()
                {
                    if (m_fileHandle != NULL)
                    {
                        ::CloseHandle(m_fileHandle);
                        m_fileHandle = NULL;
                    }

                    if (m_fileHandle2 != NULL)
                    {
                        ::CloseHandle(m_fileHandle2);
                        m_fileHandle2 = NULL;
                    }
                }

                void Reset()
                {
                    m_startPoints.clear();
                    m_exNodeCheckStatus.Clear();
                }

                std::vector<std::pair<int, float>> m_startPoints;

                CountVector<unsigned char> m_exNodeCheckStatus;

                std::vector<PageBuffer<std::uint8_t>> m_pageBuffers;

                PageBuffer<std::uint8_t> m_vectorBuffer;

                HANDLE m_fileHandle;

                HANDLE m_fileHandle2;

            };


            template<typename ValueType>
            class IExtraSearcher
            {
            public:
                IExtraSearcher()
                {
                }


                virtual ~IExtraSearcher()
                {
                }


                virtual void InitWorkSpace(ExtraWorkSpace& p_space, int p_resNumHint) = 0;

                virtual void Setup(Options& p_config) = 0;

                virtual void FinishPrepare()
                {
                }

                virtual void Search(ExtraWorkSpace& p_exWorkSpace,
                    COMMON::QueryResultSet<ValueType>& p_queryResults,
                    shared_ptr<VectorIndex> p_index,
                    SearchStats& p_stats) = 0;
            };
        }
    }
}