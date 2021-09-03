#include "inc/SSDServing/SelectHead_BKT/AnalyzeTree.h"

namespace SPTAG {
    namespace SSDServing {
        namespace SelectHead_BKT {

            void CalcLeafSize(int p_nodeID,
                const std::shared_ptr<COMMON::BKTree> p_tree,
                std::unordered_map<int, int>& p_counter) {

                SPTAG::COMMON::BKTNodeUpdate& node = (*p_tree)[p_nodeID];

                p_counter[node.centerid] = 1;

                if (node.firstChild < 0)
                {
                    return;
                }

                for (SPTAG::SizeType i = node.firstChild; i > 0; i = (*p_tree)[i].sibling)
                {
                    CalcLeafSize(i, p_tree, p_counter);
                    p_counter[node.centerid] += p_counter[(*p_tree)[i].centerid];
                }
            }

            void DfsAnalyze(int p_nodeID,
                const std::shared_ptr<COMMON::BKTree> p_tree,
                std::shared_ptr<VectorSet> p_vectorSet,
                const Options& p_opts,
                int p_height,
                std::vector<BKTNodeInfo>& p_nodeInfos) {
            
                const auto& node = (*p_tree)[p_nodeID];

                // Leaf.
                if (node.firstChild < 0)
                {
                    p_nodeInfos[p_nodeID].leafSize = 1;
                    p_nodeInfos[p_nodeID].minDepth = 0;
                    p_nodeInfos[p_nodeID].maxDepth = 0;

                    return;
                }

                p_nodeInfos[p_nodeID].leafSize = 0;

                int& minDepth = p_nodeInfos[p_nodeID].minDepth;
                int& maxDepth = p_nodeInfos[p_nodeID].maxDepth;

                int sinlgeCount = 0;
                int childrenCount = 0;
                for (int nodeId = node.firstChild; nodeId > 0; nodeId = (*p_tree)[nodeId].sibling)
                {
                    childrenCount++;
                    DfsAnalyze(nodeId, p_tree, p_vectorSet, p_opts, p_height + 1, p_nodeInfos);
                    if (minDepth > p_nodeInfos[nodeId].minDepth)
                    {
                        minDepth = p_nodeInfos[nodeId].minDepth;
                    }

                    if (maxDepth < p_nodeInfos[nodeId].maxDepth)
                    {
                        maxDepth = p_nodeInfos[nodeId].maxDepth;
                    }

                    if (p_nodeInfos[nodeId].maxDepth == 1 && p_nodeInfos[nodeId].leafSize == 1)
                    {
                        ++sinlgeCount;
                    }

                    p_nodeInfos[p_nodeID].leafSize += p_nodeInfos[nodeId].leafSize;
                }

                ++minDepth;
                ++maxDepth;

                if (p_height > 5 || sinlgeCount == 0)
                {
                    return;
                }

                LOG(Helper::LogLevel::LL_Info,
                    "CheckNode: %8d, Height: %3d, MinDepth: %3d, MaxDepth: %3d, Children: %3d, Single: %3d\n",
                    p_nodeID,
                    p_height,
                    minDepth,
                    maxDepth,
                    childrenCount,
                    sinlgeCount);

                for (int nodeId = node.firstChild; nodeId > 0; nodeId = (*p_tree)[nodeId].sibling)
                {
                    int ChildrenCount = 0;
                    for (SizeType i = (*p_tree)[i].firstChild; i > 0; i = (*p_tree)[i].sibling)
                    {
                        ChildrenCount++;
                    }
                    LOG(Helper::LogLevel::LL_Info,
                        "    ChildNode: %8d, MinDepth: %3d, MaxDepth: %3d, ChildrenCount: %3d, LeafCount: %3d\n",
                        nodeId,
                        p_nodeInfos[nodeId].minDepth,
                        p_nodeInfos[nodeId].maxDepth,
                        ChildrenCount,
                        p_nodeInfos[nodeId].leafSize);
                }
            }
        }
    }
}