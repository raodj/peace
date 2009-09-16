#ifndef MST_CLUSTER_CPP
#define MST_CLUSTER_CPP

//---------------------------------------------------------------------------
//
// Copyright (c) Miami University, Oxford, OHIO.
// All rights reserved.
//
// Miami University (MU) makes no representations or warranties about
// the suitability of the software, either express or implied,
// including but not limited to the implied warranties of
// merchantability, fitness for a particular purpose, or
// non-infringement.  MU shall not be liable for any damages suffered
// by licensee as a result of using, result of using, modifying or
// distributing this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of this license.
//
// Authors: Dhananjai M. Rao       raodm@muohio.edu
//
//---------------------------------------------------------------------------

#include "MSTCluster.h"
#include "Utilities.h"
#include "HashMap.h"
#include <cmath>

// A type-def for HashMap to make code cleaner.
typedef HashMap<int, MSTCluster*> ClusterMap;

MSTCluster::MSTCluster(MSTCluster* owner) : parent(owner) {
    // Nothing else to be done for now.
}

MSTCluster::~MSTCluster() {
    // Delete all the sub-clusters.
    for(size_t i = 0; (i < clusterList.size()); i++) {
        delete clusterList[i];
    }
    // Remove all entries from the list.
    clusterList.clear();
    // Remove all nodes in the members list.
    members.clear();
}

void
MSTCluster::add(const MSTNode& node) {
    members.push_back(node);
}

double
MSTCluster::makeClusters(NodeList& nodeList, const double percentile,
                         const int analysisCount, const ESTAnalyzer* analyzer) {
    // Compute the threshold based on the percentile value provided.
    const double threshold = calculateThreshold(nodeList.size(), percentile,
                                                analysisCount, analyzer);
    // Now create a hash map to track cluster for a given node
    ClusterMap clusterMap;
    // Now extract nodes from the nodeList and add it to appropriate
    // sub clusters.
    for(NodeList::const_iterator node = nodeList.begin();
        (node != nodeList.end()); node++) {
        if (analyzer->compareMetrics(threshold, (*node).getMetric())) {
            // Possibly this goes in a different cluster?
            continue;
        }
        // First determine if the parent of this node already has a
        // cluster associated with it.  If so, add this node to the
        // parent's cluster.
        MSTCluster *subCluster = clusterMap[(*node).parentIdx];
        if (subCluster == NULL) {
            // No entry for the parent. Time to create a new cluster...
            subCluster = new MSTCluster(this);
            // Add cluster to temporary look up hash_map. 
            clusterMap[(*node).parentIdx] = subCluster;
            // Add newly created cluster to permanent list of sub-clusters
            clusterList.push_back(subCluster);
            // Add the parent node as well to the cluster as it would
            // not have been added yet. To locate the parent node we
            // need to traverse backwards in the node list from the
            // current point.
            NodeList::const_iterator parent = node;
            while (parent != nodeList.begin()) {
                if (parent->estIdx == node->parentIdx) {
                    subCluster->add(*parent);
                    break;
                }
                // Check previous node next.
                parent--;
            }
        }
        ASSERT ( subCluster != NULL );
        subCluster->add(*node);
        // Make a temporary note of the cluster into which the node
        // has been added.
        clusterMap[(*node).estIdx] = subCluster;
    }
    // At the root node, place all un-clustered (or outliner)
    // nodes into a singleton cluster
    if (parent != NULL) {
        // This not the root node. Make the code a bit more
        // streamlined by returning righ away.
        return threshold;
    }

    // This is code path only for the root cluster.
    for(NodeList::iterator node = nodeList.begin();
        (node != nodeList.end()); node++) {
        if (clusterMap[(*node).estIdx] == NULL) {
            // This is a outlier node. Place it into a new cluster of
            // its own
            MSTCluster *subCluster = new MSTCluster(this);
            // Place the node in the new subcluster
            subCluster->add(*node);
            // Add newly created cluster to permanent list of sub-clusters
            clusterList.push_back(subCluster);
        }
    }
    // Remove all the nodes from the list of nodes.
    nodeList.clear();
    // Return threshold value back to the caller
    return threshold;
}

double
MSTCluster::calculateThreshold(const int nodeCount,
                               const double percentile,
                               const int analysisCount,
                               const ESTAnalyzer* analyzer) const {
    /*double totalSim    = 0;
    double totalSimSqr = 0;
    // Iterate over the set of nodes and compute total values to
    // determine mean and standard deviation.
    for(NodeList::const_iterator curr = nodeList.begin();
        (curr != nodeList.end()); curr++) {
        const float sim = (*curr).getMetric();
        totalSim    += sim;
        totalSimSqr += (sim * sim);
    }
    const double mean  = totalSim / (nodeList.size() - 1);
    const double stDev = sqrt((totalSimSqr / (nodeList.size() - 1))
                              - (mean * mean));
    // Compute the threshold based on percentile value provided.
    return mean + (stDev * percentile);*/

    if (!analysisCount) {
        // We are not using the TV heuristic
        // Return a threshold based on the analyzer in use
        if (!analyzer->getName().compare("baton")) {
            return 15;
        } else {
            return 40;
        }
    } else {
        // Threshold is based on ratio of full analysis (i.e. D2 runs)
        // to total number of comparisons made (EST count choose 2).
        // The idea is that this ratio will provide a measure of the overall
        // similarity of ESTs in the data set, and we can adjust our threshold
        // to compensate.  When the ratio is relatively high, we want the
        // threshold to be lower (to better differentiate the ESTs) and when
        // the ratio is relatively low, we want the threshold to be higher
        // (to cut down on singletons and give us clusters of meaningful size).

        int totalCmps = nodeCount * (nodeCount - 1) * 0.5;
        double ratio = ((double) analysisCount) / ((double) totalCmps);
        int threshold = (.08/ratio); // magic number, worked for current data
        // Put some hard limits on the threshold
        if (threshold < 40) threshold = 40;
        if (threshold > 130) threshold = 130;
        printf("Threshold: %d\n", threshold);
        return threshold;
    }
}

void
MSTCluster::makeMergedClusters(const int size, int* parent, bool* root) {
    MSTCluster *subCluster;
    ClusterMap clusterMap;
    for (int estIdx = 0; estIdx < size; estIdx++) {
        int rootIdx = estIdx;
        // replace this with actual find operation (w/ compression)
        // or move this method or something
        while (!root[rootIdx]) {
            rootIdx = parent[rootIdx];
        }
        ASSERT(root[rootIdx]);
        // This node has its own cluster
        if (clusterMap[rootIdx] == NULL) {
            // Create new cluster
            subCluster = new MSTCluster(this);
            clusterList.push_back(subCluster);
            clusterMap[rootIdx] = subCluster;
            // Add node to the cluster
            subCluster->add(MSTNode(-1, rootIdx, 0, 0));
        }
        if (rootIdx != estIdx) {
            // Need to add this EST to root's cluster as well
            subCluster = clusterMap[rootIdx];
            subCluster->add(MSTNode(-1, estIdx, 0, 0));
        }
    }
}

void
MSTCluster::printClusterTree(std::ostream& os,
                             const std::string& prefix) const {
    static int clusterID = 0;
    if (parent == NULL) {
        clusterID = 0;
    } else {
        clusterID++;
    }
    // First print information about this cluster itself.
    os << "Cluster #" << clusterID << " (#sub-clusters="
       << clusterList.size() << ", #members=" << members.size() << ")\n";
    os << prefix      << "  |\n";
    // First print nodes for this cluster.
    for(size_t i = 0; (i < members.size()); i++) {
        os << prefix << "  +---" << members[i].getESTInfo()
           << " (" << members[i].estIdx << ")\n";
    }
    // Now print each sub-cluster
    const int IterSize = clusterList.size() - 1;
    for(int i = 0; (i < IterSize); i++) {
        os << prefix << "  +---";
        clusterList[i]->printClusterTree(os, prefix + "  |   ");
    }
    if (IterSize >= 0) {
        os << prefix << "  +---";
        clusterList[IterSize]->printClusterTree(os, prefix + "      ");
    }
}

std::ostream&
operator<<(std::ostream& os, const MSTCluster& cluster) {
    static int clusterID = 0;
    clusterID = (cluster.parent != NULL) ? (clusterID + 1) : 0;
    // Print information abou this cluster.
    os << "Cluster #" << clusterID << "\n";
    for(size_t i = 0; (i < cluster.members.size()); i++) {
        os << cluster.members[i].getESTInfo() << "\n";
    }
    // Get sub-clusters to print their own information.
    for(size_t i = 0; (i < cluster.clusterList.size()); i++) {
        os << *cluster.clusterList[i] << std::endl;
    }
    return os;
}

#endif
