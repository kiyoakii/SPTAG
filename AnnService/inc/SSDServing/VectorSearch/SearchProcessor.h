#pragma once

#include "inc/SSDServing/VectorSearch/SearchStats.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/SSDServing/VectorSearch/Options.h"

namespace SPTAG {
    namespace SSDServing {
        namespace VectorSearch {
            template<typename ValueType>
            class SearchProcessor
            {
            public:
                SearchProcessor()
                {
                }


                virtual ~SearchProcessor()
                {
                }


                virtual void Setup(Options& p_config) = 0;

                virtual void Search(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats) = 0;

                virtual void SearchAsync(COMMON::QueryResultSet<ValueType>& p_queryResults, SearchStats& p_stats, std::function<void()> p_callback) = 0;

                virtual void SetHint(int p_threadNum, int p_resultNum, bool p_asyncCall, const Options& p_opts)
                {
                }
            };
        }
    }
}