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
MSTCluster::makeClusters(NodeList& nodeList, const double percentile) {
    // Compute the threshold based on the percentile value provided.
    const double threshold = calculateThreshold(nodeList, percentile);
    // Now create a hash map to track cluster for a given node
    ClusterMap clusterMap;
    // Now extract nodes from the nodeList and add it to appropriate
    // sub clusters.
    for(NodeList::const_iterator node = nodeList.begin();
        (node != nodeList.end()); node++) {
        if ((*node).getMetric() >= threshold) {
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
MSTCluster::calculateThreshold(const NodeList& UNREFERENCED_PARAMETER(nodeList),
                               const double UNREFERENCED_PARAMETER(percentile)) const {
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

    return 130;
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
