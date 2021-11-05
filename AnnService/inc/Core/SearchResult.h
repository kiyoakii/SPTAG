// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#ifndef _SPTAG_SEARCHRESULT_H_
#define _SPTAG_SEARCHRESULT_H_

#include "CommonDataStructure.h"

namespace SPTAG
{
    struct BasicResult
    {
        SizeType VID;
        float Dist;
        ByteArray Meta;

        SizeType fatherVID;

        BasicResult() : VID(-1), Dist(MaxDist), fatherVID(-1) {}

        BasicResult(SizeType p_vid, float p_dist) : VID(p_vid), Dist(p_dist), fatherVID(-1) {}
        
        BasicResult(SizeType p_vid, float p_dist, SizeType f_vid) : VID(p_vid), Dist(p_dist), fatherVID(f_vid) {}

        BasicResult(SizeType p_vid, float p_dist, ByteArray p_meta) : VID(p_vid), Dist(p_dist), Meta(p_meta), fatherVID(-1) {}
    };

} // namespace SPTAG

#endif // _SPTAG_SEARCHRESULT_H_
