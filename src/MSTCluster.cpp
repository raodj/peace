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
#include "arg_parser.h"
#include <cmath>

// The static variable to generate unique cluster ID values.
int MSTCluster::clusterIDSequence = 0;

// A type-def for HashMap to make code cleaner.
typedef HashMap<int, MSTCluster*> ClusterMap;

MSTCluster::MSTCluster(MSTCluster* owner) :
    parent(owner), clusterID(clusterIDSequence++) {
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
MSTCluster::makeClusters(NodeList& nodeList, const ESTAnalyzer* analyzer) {
    const float threshold = 130;
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
	std::string newPrefix = prefix;
	newPrefix += "  |   ";
        clusterList[i]->printClusterTree(os, newPrefix);
    }
    if (IterSize >= 0) {
        os << prefix << "  +---";
	std::string newPrefix = prefix;
	newPrefix += "      ";
        clusterList[IterSize]->printClusterTree(os, newPrefix);
    }
}

void
MSTCluster::guiPrintClusterTree(std::ostream& os, const char *srcFile) const {
    // This  method must be called only on the root node.
    ASSERT ( parent == NULL );
    // Buffer for string of time for reporting
    char  now[128], srcTimeStr[128] = "<none>";
    if (srcFile != NULL) {
        // Get time stamp of source file for coherence checks.
        getTimeStamp(srcFile, srcTimeStr);
    }
    // Dump meta data about how this cluster was built.
    os << "# Cluster Data\n"
       << "# Generated on: "    << getTime(now)
       << "# Generated from source file: "
       << ((srcFile != NULL) ? srcFile : "<none>") << "\n"
       << "# Source file timestamp: " << srcTimeStr
       << "# Data format: <E,estIdx,clstrId> | <C,clstrID,ParntClstrID>\n"
       << "# Command Line: " << arg_parser::get_global_args() << "\n";
    // Now print the mst cluster for GUI processing
    guiPrintTree(os);
}

void
MSTCluster::guiPrintTree(std::ostream& os) const {
    // Print information about this cluster and its parent cluster.
    os << "C," << clusterID << ","
       << ((parent != NULL) ? parent->clusterID : -1) << std::endl;
    // Print all the ESTs in this cluster.
    for(size_t i = 0; (i < members.size()); i++) {
        os << "E," << members[i].getESTIdx()
           << ","  << clusterID << "\n";
    }
    // Get sub-clusters to print their own information.
    for(size_t i = 0; (i < clusterList.size()); i++) {
        clusterList[i]->guiPrintTree(os);
    }
}

std::ostream&
operator<<(std::ostream& os, const MSTCluster& cluster) {
    // Print information abou this cluster.
    os << "Cluster #" << cluster.clusterID << "\n";
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
