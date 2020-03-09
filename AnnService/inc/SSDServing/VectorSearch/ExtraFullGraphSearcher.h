#pragma once
#include "inc/SSDServing/Common/stdafx.h"
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/SSDServing/VectorSearch/SearchStats.h"

#include <fstream>
#include <map>
#include <memory>
#include <vector>
#include <future>

#include <fileapi.h>

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {
            void ErrorExit()
            {
                // Retrieve the system error message for the last-error code

                LPVOID lpMsgBuf;
                DWORD dw = GetLastError();

                FormatMessage(
                    FORMAT_MESSAGE_ALLOCATE_BUFFER |
                    FORMAT_MESSAGE_FROM_SYSTEM |
                    FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL,
                    dw,
                    0,
                    (LPTSTR)&lpMsgBuf,
                    0, NULL);

                // Display the error message and exit the process

                std::fprintf(stderr, "Failed with: %s\n", (char*)lpMsgBuf);

                LocalFree(lpMsgBuf);
                ExitProcess(dw);
            }

            template <typename ValueType>
            class ExtraFullGraphSearcher : public IExtraSearcher<ValueType>
            {
            public:
                ExtraFullGraphSearcher(const std::string& p_extraFullGraphFile)
                {
                    m_extraFullGraphFile = p_extraFullGraphFile;
                    LoadingHeadInfo(p_extraFullGraphFile);
                }

                virtual ~ExtraFullGraphSearcher()
                {
                }

                virtual void InitWorkSpace(ExtraWorkSpace& p_space, int p_resNumHint)
                {
                    p_space.m_startPoints.clear();
                    p_space.m_startPoints.reserve(p_resNumHint);

                    p_space.m_exNodeCheckStatus.Init(m_totalDocumentCount);

                    p_space.m_pageBuffers.resize(128);
                    for (auto& pageBuffer : p_space.m_pageBuffers)
                    {
                        pageBuffer.ReservePageBuffer(c_pageSize * 2);
                    }

                    p_space.m_fileHandle = ::CreateFileA(m_extraFullGraphFile.c_str(),
                        GENERIC_READ,
                        FILE_SHARE_READ,
                        NULL,
                        OPEN_EXISTING,
                        FILE_FLAG_NO_BUFFERING,
                        NULL);

                    if (p_space.m_fileHandle == INVALID_HANDLE_VALUE)
                    {
                        fprintf(stderr, "Failed to create file handle: %s\n", m_extraFullGraphFile.c_str());
                        exit(1);
                    }

                    p_space.m_fileHandle2 = ::CreateFileA(m_extraFullGraphFile.c_str(),
                        GENERIC_READ,
                        FILE_SHARE_READ,
                        NULL,
                        OPEN_EXISTING,
                        FILE_FLAG_NO_BUFFERING,
                        NULL);

                    if (p_space.m_fileHandle2 == INVALID_HANDLE_VALUE)
                    {
                        fprintf(stderr, "Failed to create file handle: %s\n", m_extraFullGraphFile.c_str());
                        exit(1);
                    }
                }

                virtual void Setup(Options& p_config)
                {
                    m_maxCheck = atoi(p_config.m_extraMaxCheck.c_str());

                    m_parallelLoadPercentage = atoi(p_config.m_parallelLoadPercentage.c_str());
                    m_parallelLoadPercentage = std::min<int>(m_parallelLoadPercentage, 100);
                    m_parallelLoadPercentage = std::max<int>(m_parallelLoadPercentage, 0);
                }

                virtual void Search(ExtraWorkSpace& p_exWorkSpace,
                    COMMON::QueryResultSet<ValueType>& p_queryResults,
                    shared_ptr<VectorIndex> p_index,
                    SearchStats& p_stats)
                {
                    int dedupeCount = 0;
                    if (p_exWorkSpace.m_pageBuffers.size() < p_exWorkSpace.m_startPoints.size())
                    {
                        p_exWorkSpace.m_pageBuffers.resize(p_exWorkSpace.m_startPoints.size());
                    }

                    bool finish = false;
                    int curCheck = 0;
                    size_t halfP = static_cast<size_t>(p_exWorkSpace.m_startPoints.size() * (100 - m_parallelLoadPercentage) / 100);

                    std::atomic_int32_t diskRead = 0;

                    auto loadHalf = [&]()
                    {
                        for (size_t i = halfP; i < p_exWorkSpace.m_startPoints.size(); ++i)
                        {
                            const auto& vi = p_exWorkSpace.m_startPoints[i];
                            if (m_listInfos[vi.first].listEleCount == 0)
                            {
                                continue;
                            }

                            diskRead += m_listInfos[vi.first].listPageCount;
                            ReadList(p_exWorkSpace.m_fileHandle2, p_exWorkSpace.m_pageBuffers[i], m_listInfos[vi.first]);
                        }
                    };

                    std::future<void> loadHalfThread = std::async(std::launch::async, loadHalf);

                    // float currMaxDist = std::sqrtf(p_queryResults.GetResult(0).Score);
                    auto traverseLists = [&](size_t p_si, size_t p_ei, bool p_needRead)
                    {
                        for (size_t si = p_si; si < p_ei && !finish; ++si)
                        {
                            const auto& vi = p_exWorkSpace.m_startPoints[si];
                            if (m_listInfos[vi.first].listEleCount == 0)
                            {
                                continue;
                            }

                            // float headQueryDist = std::sqrtf(vi.second);

                            if (p_needRead)
                            {
                                diskRead += m_listInfos[vi.first].listPageCount;
                                ReadList(p_exWorkSpace.m_fileHandle, p_exWorkSpace.m_pageBuffers[si], m_listInfos[vi.first]);
                            }

                            for (int i = 0; i < m_listInfos[vi.first].listEleCount; ++i)
                            {
                                std::uint8_t* vectorInfo = p_exWorkSpace.m_pageBuffers[si].GetBuffer() + m_listInfos[vi.first].pageOffset + i * m_vectorInfoSize;
                                int vectorID = *(reinterpret_cast<int*>(vectorInfo));
                                vectorInfo += sizeof(int);
                                // float distToHead = *(reinterpret_cast<float*>(vectorInfo));
                                // vectorInfo += sizeof(float);

                                if (p_exWorkSpace.m_exNodeCheckStatus.CheckAndSet(vectorID))
                                {
                                    ++dedupeCount;
                                    continue;
                                }

                                /*
                                if (distToHead < (headQueryDist - currMaxDist) * 0.8 || distToHead > (headQueryDist + currMaxDist) * 1.2)
                                {
                                    ++p_stats.m_distCheckFilterCount;
                                    continue;
                                }
                                */

                                ++curCheck;

                                auto distance2leaf = p_index->ComputeDistance(p_queryResults.GetTarget(),
                                    vectorInfo);

                                if (p_queryResults.AddPoint(vectorID, distance2leaf))
                                {
                                    //  currMaxDist = std::sqrtf(p_queryResults.GetResult(0).Score);
                                }

                                if (curCheck >= m_maxCheck)
                                {
                                    finish = true;
                                    break;
                                }
                            }
                        }
                    };

                    traverseLists(0, halfP, true);
                    loadHalfThread.wait();
                    traverseLists(halfP, p_exWorkSpace.m_startPoints.size(), false);

                    p_stats.m_exCheck = curCheck;
                    p_stats.m_diskAccessCount = diskRead;

                    p_queryResults.SortResult();
                }


            private:
                struct ListInfo
                {
                    int pageNum = 0;

                    std::uint16_t pageOffset = 0;

                    int listEleCount = 0;

                    std::uint16_t listPageCount = 0;
                };

                void ReadList(HANDLE& p_fileHandle, PageBuffer<std::uint8_t>& p_buffer, const ListInfo& p_listInfo)
                {
                    size_t totalBytes = static_cast<size_t>(p_listInfo.listPageCount)* c_pageSize;

                    p_buffer.ReservePageBuffer(totalBytes);

                    DWORD bytesRead = 0;

                    LARGE_INTEGER li;
                    li.QuadPart = (m_listPageOffset + p_listInfo.pageNum) * c_pageSize;

                    if (!::SetFilePointerEx(p_fileHandle, li, NULL, FILE_BEGIN))
                    {
                        ErrorExit();
                    }

                    if (!::ReadFile(p_fileHandle,
                        p_buffer.GetBuffer(),
                        static_cast<DWORD>(totalBytes),
                        &bytesRead,
                        NULL))
                    {
                        ErrorExit();
                    }
                }

            private:
                void LoadingHeadInfo(const std::string& p_file)
                {
                    std::ifstream input(p_file, std::ios::binary);
                    if (!input.is_open())
                    {
                        fprintf(stderr, "Failed to open file: %s\n", p_file.c_str());
                        exit(1);
                    }

                    input.read(reinterpret_cast<char*>(&m_listCount), sizeof(m_listCount));
                    input.read(reinterpret_cast<char*>(&m_totalDocumentCount), sizeof(m_totalDocumentCount));
                    input.read(reinterpret_cast<char*>(&m_iDataDimension), sizeof(m_iDataDimension));
                    input.read(reinterpret_cast<char*>(&m_listPageOffset), sizeof(m_listPageOffset));

                    m_listInfos.reset(new ListInfo[m_listCount]);

                    size_t totalListElementCount = 0;

                    std::map<int, int> pageCountDist;

                    int biglistCount = 0;
                    int biglistElementCount = 0;
                    for (int i = 0; i < m_listCount; ++i)
                    {
                        input.read(reinterpret_cast<char*>(&(m_listInfos[i].pageNum)), sizeof(m_listInfos[i].pageNum));
                        input.read(reinterpret_cast<char*>(&(m_listInfos[i].pageOffset)), sizeof(m_listInfos[i].pageOffset));
                        input.read(reinterpret_cast<char*>(&(m_listInfos[i].listEleCount)), sizeof(m_listInfos[i].listEleCount));
                        input.read(reinterpret_cast<char*>(&(m_listInfos[i].listPageCount)), sizeof(m_listInfos[i].listPageCount));

                        totalListElementCount += m_listInfos[i].listEleCount;
                        int pageCount = m_listInfos[i].listPageCount;

                        if (pageCount > 1)
                        {
                            ++biglistCount;
                            biglistElementCount += m_listInfos[i].listEleCount;
                        }

                        if (pageCountDist.count(pageCount) == 0)
                        {
                            pageCountDist[pageCount] = 1;
                        }
                        else
                        {
                            pageCountDist[pageCount] += 1;
                        }
                    }

                    input.close();

                    fprintf(stderr,
                        "Finish reading header info, list count %d, total doc count %d, dimension %d, list page offset %d.\n",
                        m_listCount,
                        m_totalDocumentCount,
                        m_iDataDimension,
                        m_listPageOffset);


                    fprintf(stderr,
                        "Big page (>4K): list count %d, total element count %d.\n",
                        biglistCount,
                        biglistElementCount);

                    fprintf(stderr, "Total Element Count: %llu\n", totalListElementCount);

                    for (auto& ele : pageCountDist)
                    {
                        fprintf(stderr, "Page Count Dist: %d %d\n", ele.first, ele.second);
                    }

                    m_vectorInfoSize = m_iDataDimension * sizeof(ValueType) + sizeof(int);
                }


            private:

                std::string m_extraFullGraphFile;

                const static std::uint64_t c_pageSize = 4096;

                int m_parallelLoadPercentage;

                int m_maxCheck;

                int m_iDataDimension;

                int m_listCount;

                int m_totalDocumentCount;

                int m_listPageOffset;

                size_t m_vectorInfoSize;

                std::unique_ptr<ListInfo[]> m_listInfos;
            };
        }
    }
}
