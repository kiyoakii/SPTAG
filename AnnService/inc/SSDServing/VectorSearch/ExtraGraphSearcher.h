#pragma once
#include <memory>
#include <vector>
#include "inc/SSDServing/VectorSearch/IExtraSearcher.h"
#include "inc/Helper/DynamicNeighbors.h"

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {
            template<typename ValueType>
            class ExtraGraphSearcher : public IExtraSearcher<ValueType>
            {
            public:
                ExtraGraphSearcher(const std::string p_extraVectorSetFile,
                    const std::string& p_neighborsFile)
                    : m_maxCheck(1024)
                {
                    fprintf(stderr, "Loading extra graph vector set...\n");
                    m_extraVectorSet.reset(new BasicVectorSet(p_extraVectorSetFile.c_str(), GetEnumValueType<ValueType>()));

                    fprintf(stderr, "Loading extra graph set...\n");
                    m_dynamicNeighborsSet.reset(new Helper::DynamicNeighborsSet(p_neighborsFile.c_str()));

                    fprintf(stderr, "Finish loading extra graph set.\n");
                }

                virtual ~ExtraGraphSearcher()
                {
                }


                virtual void InitWorkSpace(ExtraWorkSpace* p_space, int p_resNumHint)
                {
                    p_space->m_deduper.Clear();
                }


                virtual void Search(ExtraWorkSpace* p_exWorkSpace,
                    COMMON::QueryResultSet<ValueType>& p_queryResults,
                    shared_ptr<VectorIndex> p_index,
                    SearchStats& p_stats)
                {
                    InitWorkSpace(p_exWorkSpace, p_exWorkSpace->m_postingIDs.size());

                    bool finish = false;
                    int curCheck = 0;

                    for (const auto& vi : p_exWorkSpace->m_postingIDs)
                    {
                        const auto& nl = (*m_dynamicNeighborsSet)[vi];

                        p_stats.m_totalListElementsCount += static_cast<int>(nl.Size());

                        for (int i = 0; i < nl.Size(); ++i)
                        {
                            int v = nl[i];
                            if (p_exWorkSpace->m_deduper.CheckAndSet(v))
                            {
                                continue;
                            }

                            ++curCheck;

                            auto distance2leaf = p_index->ComputeDistance(p_queryResults.GetTarget(),
                                m_extraVectorSet->GetVector(v));

                            p_queryResults.AddPoint(v, distance2leaf);
                            if (curCheck >= m_maxCheck)
                            {
                                finish = true;
                                break;
                            }
                        }

                        if (finish)
                        {
                            break;
                        }
                    }

                    p_stats.m_exCheck = curCheck;
                    p_queryResults.SortResult();
                }

                virtual void Setup(Options& p_config)
                {
                    m_maxCheck = atoi(p_config.m_extraMaxCheck.c_str());
                }

            private:
                std::unique_ptr<BasicVectorSet> m_extraVectorSet;

                std::unique_ptr<Helper::DynamicNeighborsSet> m_dynamicNeighborsSet;

                int m_maxCheck;
            };
        }
    }
}
