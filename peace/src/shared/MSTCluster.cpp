#ifndef MST_CLUSTER_CPP
#define MST_CLUSTER_CPP

//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "MSTCluster.h"
#include "Utilities.h"
#include "HashMap.h"
#include "ESTList.h"
#include <cmath>
#include <sstream>

// The static variable to generate unique cluster ID values.
int MSTCluster::clusterIDSequence = 0;

// A type-def for HashMap to make code cleaner.
typedef HashMap<int, MSTCluster*> ClusterMap;

// Global list to maintain reference to all the MSTClusters created.
ClusterList MSTCluster::globalClusterList;

MSTCluster::MSTCluster(MSTCluster* owner, const std::string& clsName) :
    parent(owner), clusterID(clusterIDSequence++), name(clsName) {
    // Save reference to this cluster entry in the global list
    ASSERT ( clusterID == (int) globalClusterList.size() );
    globalClusterList.push_back(this);
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
    // Clear out entry in the global EST list
    globalClusterList[clusterID] = NULL;
}

void
MSTCluster::add(const MSTNode& node) {
    members.push_back(node);
}

void
MSTCluster::add(const int clusterID, const MSTNode& node) {
    // Validate context and parameters
    ASSERT ((clusterID >= 0) && (clusterID < (int) globalClusterList.size()));
    ASSERT (globalClusterList[clusterID] != NULL);
    ASSERT (globalClusterList[clusterID]->clusterID == clusterID);
    // Add node to the specified cluster.
    globalClusterList[clusterID]->members.push_back(node);
}

void
MSTCluster::add(MSTCluster *child) {
    ASSERT ( child->parent == this );
    // Add newly created cluster to list of sub-clusters
    clusterList.push_back(child);
}

void
MSTCluster::makeClusters(const NodeList& nodeList, const ESTAnalyzer* analyzer,
                         float threshold) {
    // Check if the threshold is less than 0 (treat all negatives as -1)
    if (threshold < 0) {
        threshold = calculateThreshold(nodeList);
    }
    // Create a hash map to track cluster for a given node
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
            // Add to list of sub-clusters
            add(subCluster);
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
        // streamlined by returning right away.
        return;
    }

    // This is code path only for the root cluster.
    for(NodeList::const_iterator node = nodeList.begin();
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
}

float
MSTCluster::calculateThreshold(const NodeList& nodeList) const {
    double totalSim = 0;
    double totalSimSqr = 0;
    // Iterate over the set of nodes and compute total values
    // to determine mean and standard deviation.
    for (NodeList::const_iterator curr = nodeList.begin();
         (curr != nodeList.end()); curr++) {
        const float sim = (*curr).getMetric();
        totalSim    += sim;
        totalSimSqr += (sim*sim);
    }
    const double mean  = totalSim / (nodeList.size()-1);
    const double stDev = sqrt((totalSimSqr / (nodeList.size()-1))
                              - (mean*mean));
    // Return the threshold
    return (float) (mean + stDev);
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
MSTCluster::printClusterTree(const ESTList& estList,
                             std::ostream& os,
                             const std::string& prefix) const {
    // First print information about this cluster itself.
    os << "Cluster #" << clusterID;
    if (name != "") {
        os << " [" << name << "]";
    }
    os << " (#sub-clusters="
       << clusterList.size() << ", #members=" << members.size() << ")\n";
    os << prefix      << "  |\n";
    // First print nodes for this cluster.
    for(size_t i = 0; (i < members.size()); i++) {
        os << prefix << "  +---" << getESTInfo(members[i].estIdx, estList)
           << " (" << members[i].estIdx << ")\n";
    }
    // Now print each sub-cluster
    const int IterSize = clusterList.size() - 1;
    for(int i = 0; (i < IterSize); i++) {
        os << prefix << "  +---";
        std::string newPrefix = prefix;
        newPrefix += "  |   ";
        clusterList[i]->printClusterTree(estList, os, newPrefix);
    }
    if (IterSize >= 0) {
        os << prefix << "  +---";
        std::string newPrefix = prefix;
        newPrefix += "      ";
        clusterList[IterSize]->printClusterTree(estList, os, newPrefix);
    }
}

void
MSTCluster::guiPrintClusterTree(std::ostream& os) const {
    // This  method must be called only on the root node.
    ASSERT ( parent == NULL );
    // Now print the mst cluster for GUI processing
    guiPrintTree(os);
}

void
MSTCluster::guiPrintTree(std::ostream& os) const {
    // Print information about this cluster and its parent cluster.
    os << "C," << clusterID << ","
       << ((parent != NULL) ? parent->clusterID : -1) 
       << "," << name << std::endl;
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
    os << "Cluster #" << cluster.clusterID << " ["
       << cluster.name << "]\n";
    for(size_t i = 0; (i < cluster.members.size()); i++) {
        os << cluster.members[i] << "\n";
    }
    // Get sub-clusters to print their own information.
    for(size_t i = 0; (i < cluster.clusterList.size()); i++) {
        os << *cluster.clusterList[i] << std::endl;
    }
    return os;
}

std::string
MSTCluster::getESTInfo(const int estIdx, const ESTList& estList) const {
    std::string info = "";
    const EST *est   = estList.get(estIdx);
    if ((est != NULL) && (est->getInfo() != NULL)) {
        info = est->getInfo();
    } else {
	std::ostringstream estInfo;
	estInfo << "ESTidx #" << estIdx;
	info = estInfo.str();
    }
    return info;
}

#endif
