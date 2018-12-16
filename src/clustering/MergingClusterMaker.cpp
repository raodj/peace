#ifndef MERGING_CLUSTER_MAKER_CPP
#define MERGING_CLUSTER_MAKER_CPP

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include <algorithm>
#include "MergingClusterMaker.h"
#include "MPIHelper.h"
#include "MPIStats.h"
#include "ESTList.h"
#include "HeuristicChain.h"

// The static analyzer pointer in ClsInfo inner class
ESTAnalyzer* MergingClusterMaker::ClsInfo::analyzer = NULL;

MergingClusterMaker::MergingClusterMaker(ESTAnalyzer *analyzer)
    : MSTClusterMaker("merging", analyzer) {
    // Initialize other instance variables & command-line arguments to
    // custom default values.
    mergeThresh = 0.01;  // 1%
    mergeStride = 0.1;   // Check max 10% of reads in each cluster
    mergeReps   = -1;    // Keep merging until no change

    // The default heuristic to be removed prior to comparisons
    rmHeuristics.push_back("primes");
}

MergingClusterMaker::~MergingClusterMaker() {
    // Nothing to be done for now. Base class cleans up some of the
    // dynamic memory-based instance variables.
}

// Add additional command-line options for this cluster maker.
void
MergingClusterMaker::addCommandLineArguments(ArgParser& argParser) {
    // Let base class do its thing first
    MSTClusterMaker::addCommandLineArguments(argParser);
    // Now define our local command line arguments
    const ArgParser::ArgRecord LocalArgs[] = {
        {"--cls-merge-thresh", "Cluster size (in % of reads) to be merged",
         &mergeThresh, ArgParser::FLOAT},
        {"--cls-merge-stride", "% of reads to compare to merge clusters",
         &mergeStride, ArgParser::FLOAT},    
        {"--cls-merge-reps", "Maximum tries to merge clusters",
         &mergeReps, ArgParser::BOOLEAN},
        {"--cls-merge-rm-heur", "List of heuristics to be ignored",
         &rmHeuristics, ArgParser::STRING_LIST},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add our command line arguments to the parser.
    argParser.addValidArguments(LocalArgs);
}

int
MergingClusterMaker::makeClusters() {
    // First let the base class create the initial set of clusters to
    // be merged.
    int errCode = MSTClusterMaker::makeClusters();
    if (errCode != NO_ERROR) {
        return errCode;  // Error occurred during initial clustering 
    }
    // Setup EST analyzer reference in our inner class
    ClsInfo::analyzer = analyzer;
    // NOTE: The clustering information is available only on the
    // manager.  So this information needs to be broadcasted out to
    // all the worker processes so they all have consistent
    // information.
    
    // First broadcast the clustering threshold so that all the
    // processes have the same threshold to work with.
    MPI_BCAST(&clsThreshold, 1, MPI_TYPE_FLOAT, MANAGER_RANK);

    // Compute the reference ESTs (list is built only on the manager).
    IntIntList refESTs = getRefESTs();
    // Broadcast list to all work processes (in 2 steps)
    int refListSize = refESTs.size();
    MPI_BCAST(&refListSize, 1, MPI_TYPE_INT, MANAGER_RANK);
    refESTs.resize(refListSize);  // Ensure sufficient space
    MPI_BCAST(&refESTs[0], refListSize, MPI_TYPE_2INT, MANAGER_RANK);

    // Now the manager and worker processes have the necessary
    // information to collectively start merging clusters.  Next,
    // remove heuristics to be ignored from heuristic chain.
    for (size_t i = 0; (i < rmHeuristics.size()); i++) {
        analyzer->getHeuristicChain()->removeHeuristic(rmHeuristics[i]);
    }

    // Figure out the sub-set of reads that this process should work
    // on when checking for matching clusters
    int startIndex = 0, endIndex = 0;
    ::getLocallyOwnedESTidx(refListSize, startIndex, endIndex);

    // Now, do different operations on manager and worker processes --
    // workers help to find closely matching clusters.  Manager does
    // all the work of identifying and merging clusters.
    if (MPI_GET_RANK() == MANAGER_RANK) {
        errCode = mergeManager(refESTs, startIndex, endIndex);
    } else {
        errCode = mergeWorker(refESTs, startIndex, endIndex);
    }
    // Everything went well.
    return errCode;
}

int
MergingClusterMaker::mergeManager(const IntIntList& refESTs,
                                  const int startIndex, const int endIndex) {
    // Get the size of all the clusters (sorted in ascending order) to
    // determine which ones need to be merged.
    IntIntList clsSizes = getClusterSizes();
    // Find out the threshold (in number of reads) below which
    // clusters should be collapsed.
    const int sizeThresh = std::max<int>(1, estList->size() * mergeThresh);
    // Now check each cluster and merge it if needed.
    for (size_t clsIdx = 0; (clsIdx < clsSizes.size()); clsIdx++) {
        if (clsSizes[clsIdx].second <= sizeThresh) {
            // This cluster size is below the size threshold. Check to
            // see if it can be merged.
            const int clsId       = clsSizes[clsIdx].first;
            const int targetClsId = findBestCluster(clsId, refESTs,
                                                    startIndex, endIndex);
            if (targetClsId != -1) {
                ASSERT(targetClsId != clsId);
                // Found a cluster into which the this cluster can be
                // merged. So do the merge.
                MSTCluster* target = MSTCluster::getCluster(targetClsId);
                ASSERT( target != NULL );
                target->mergeCluster(clsId);
                // Let workers know that a cluster has been merged and
                // this cluster shoudld be ignored in future
                // comparisons.
                mergedClusterIDs.insert(clsId);
                IntInt info(clsId, -1);
                MPI_BCAST(&info, 1, MPI_TYPE_2INT, MANAGER_RANK);
                std::cerr << "Merged cluster #" << clsId
                          << " into #" << targetClsId << std::endl;
            } else {
                std::cerr << "Did not find a matching cluster for #"
                          << clsId << std::endl;
            }
        }
    }
    // Finally send sentinel values to workers to let them know we are
    // all done.
    IntInt info(-1, -1);
    MPI_BCAST(&info, 1, MPI_TYPE_2INT, MANAGER_RANK);
    // All done successfully.
    return NO_ERROR;
}

int
MergingClusterMaker::findBestCluster(const int clsId, const IntIntList& refESTs,
                                     const int startIndex, const int endIndex) {
    // This method should be called only the manager process.
    ASSERT(MPI_GET_RANK() == MANAGER_RANK);
    // Get the reference reads for the given cluster
    IntIntList clsRefESTs = getRefESTs(clsId);
    // Check each reference read in the cluster until we find a
    // matching cluster
    for (size_t estIdx = 0; (estIdx < clsRefESTs.size()); estIdx++) {
        // Broadcast the reference read and reference cluster ID to the
        // worker processes.
        IntInt clsEst(clsId, clsRefESTs[estIdx].first);
        MPI_BCAST(&clsEst, 1, MPI_TYPE_2INT, MANAGER_RANK);
        // Compute the local matching clusters.
        ClsInfoList clsList = findMatchingClusters(clsEst.second, clsEst.first,
                                                   refESTs, startIndex,
                                                   endIndex);
        // Collect sorted list of values from all the workers
        mergeFromWorkers(clsList);
        // Sort the list to place the the best matching cluster at the top
        std::sort(clsList.begin(), clsList.end());
        // Check if we have a good match
        if (!clsList.empty() &&
            analyzer->compareMetrics(clsList[0].metric, clsThreshold)) {
            // Yes! we found a cluster to merge into.
            return clsList[0].clusterID;
        }
    }
    // No matching cluster found.
    return -1;
}

void
MergingClusterMaker::mergeFromWorkers(ClsInfoList& clsList) {
    const int WorkerCount = MPI_GET_SIZE();
    // Get cluster information from each worker.
    for(int workerID = 1; (workerID < WorkerCount); workerID++) {
        // Wait for a message to be received from a worker and obtain
        // some status information regarding the message.
        MPI_STATUS msgInfo;
        MPI_CODE({
                const int srcRank = (strictOrder ? workerID : MPI_ANY_SOURCE);
                MPI_PROBE(srcRank, CLUSTER_INFO_LIST, msgInfo);
            });
        // OK, we have a valid repopulation request pending from some
        // worker. So read and process it.
        const int bytes = MPI_GET_COUNT(msgInfo, MPI_TYPE_CHAR);
        ASSERT(bytes >= (int) sizeof(ClsInfo));
        ASSERT(bytes % sizeof(ClsInfo) == 0);
        ClsInfoList workerList(bytes / sizeof(ClsInfo));
        MPI_RECV(&workerList[0], bytes, MPI_TYPE_CHAR, msgInfo.MPI_SOURCE,
                 CLUSTER_INFO_LIST);
        // Add valid entries from the worker into the supplied clsList
        for (size_t i = 0; (i < workerList.size()); i++) {
            if (workerList[i].clusterID != -1) {
                workerList[i].clusterSize =
                    getClusterSize(workerList[i].clusterID);
                clsList.push_back(workerList[i]);
            }
        }
    }
}

int
MergingClusterMaker::mergeWorker(const IntIntList& refESTs,
                                  const int startIndex, const int endIndex) {
    bool done = false;  // Flag to end end loop below.
    do {
        // Recieve broadcast from manager about cluster & EST IDs
        IntInt clsEstIDs;  // Sentinel values broadcast by manager
        MPI_BCAST(&clsEstIDs, 1, MPI_TYPE_2INT, MANAGER_RANK);
        // Based on cluster ID and EST ID we perform different tasks.
        done = (clsEstIDs.first == -1);
        if (!done) {
            // Based on EST ID we take two different approaches
            if (clsEstIDs.second != -1) {
                // Compute local matching clusters
                ClsInfoList clsList = findMatchingClusters(clsEstIDs.second,
                                                           clsEstIDs.first,
                                                           refESTs, startIndex,
                                                           endIndex);
                // Send them to the manager.
                MPI_SEND(&clsList[0], clsList.size() * sizeof(ClsInfo),
                         MPI_TYPE_CHAR, MANAGER_RANK, CLUSTER_INFO_LIST);
            } else {
                // A cluster has been merged. So track the merged
                // cluster so that we can ignore the reads associated
                // with it.
                mergedClusterIDs.insert(clsEstIDs.first);
            }
        }
    } while (!done);
    // Everything went well
    return NO_ERROR;
}

MergingClusterMaker::ClsInfoList
MergingClusterMaker::findMatchingClusters(const int estId,
                                          const int srcCluster,
                                          const IntIntList& refESTs,
                                          const int startIndex,
                                          const int endIndex) {
    ClsInfoList clsList;  // List populated in the loop below.
    // Set the reference EST to be used
    ASSERT(estList != NULL);
    analyzer->setReferenceEST(estList->get(estId, true));
    // Compare with the reference ESTs from other clusters to find
    // matching cluster.
    for (int i = startIndex; (i < endIndex); i++) {
        // Ignore reads in already merged clusters
        if ((refESTs[i].second == srcCluster) ||
            (mergedClusterIDs.find(refESTs[i].second) != mergedClusterIDs.end())) {
            continue;  // Cluster already merged. Ignore this read
        }
        // Get comparison metric (ignoring some heuristics)
        const float metric = analyze(estList->get(refESTs[i].first));
        ASSERT ( metric >= 0 );
        // Check if metric is below clustering threshold. 
        if (analyzer->compareMetrics(metric, clsThreshold)) {
            // Found good candidate target cluster. Track its
            // information.
            const int clsId = refESTs[i].second;
            clsList.push_back(ClsInfo(clsId, metric, getClusterSize(clsId)));
        }
    }
    // If the list is empty then we add a dummy entry to send to the
    // manager.
    if (clsList.empty() && (MPI_GET_RANK() != MANAGER_RANK)) {
        // Add an invalid entry.
        clsList.push_back(ClsInfo());
    }
    // The list of closely matching clusters.
    return clsList;
}

int
MergingClusterMaker::getClusterSize(const int clsID) const {
    const MSTCluster* const cls = MSTCluster::getCluster(clsID);
    return (cls != NULL ? cls->getMembers().size() : 0);
}

IntIntList
MergingClusterMaker::getRefESTs() const {
    IntIntList refESTs;  // vector to be populated below
    // Get the global list of clusters.
    const ClusterList& clsList = MSTCluster::getGlobalClusterList();    
    // Build list of reference reads for every cluster
    for (size_t i = 1; (i < clsList.size()); i++) {
        if (clsList[i] == NULL) {
            continue;  // not a cluster entry
        }
        // Get reference reads for this cluster via helper method
        const IntIntList subList = getRefESTs(clsList[i]->getClusterID());
        // Add sublist the main list.
        refESTs.insert(refESTs.end(), subList.begin(), subList.end());
    }
    // Return the full list of of reference reads back.
    return refESTs;
}

IntIntList
MergingClusterMaker::getRefESTs(const int clusterID) const {
    IntIntList refESTs;  // vector to be populated below
    // First, get the cluster we are processing.
    const MSTCluster* const cls = MSTCluster::getCluster(clusterID);
    if ((cls == NULL) || (cls->getMembers().empty())) {
        return refESTs;  // invalid or empty cluster ID
    }
    // Get the list of reads for this cluster
    const NodeList& members = cls->getMembers();
    size_t refIdx = 0;  // Index of reference read
    while (refIdx < members.size()) {
        // Add the current reference read.
        refESTs.push_back(IntInt(members[refIdx].getESTIdx(), clusterID));
        // Skip over reads based on the stride specified.
        refIdx += std::max<int>(1, (mergeStride * members.size()));
    }
    // The list of reference reads for this cluster.
    return refESTs;
}

IntIntList
MergingClusterMaker::getClusterSizes() const {
    IntIntList clsSizes;  // vector to be populated below.
    // Get the global list of clusters.
    const ClusterList& clsList = MSTCluster::getGlobalClusterList();
    // Process each entry in the cluster list.
    for (size_t i = 1; (i < clsList.size()); i++) {
        if (clsList[i] == NULL) {
            continue;  // not a cluster entry
        }
        // Add entry to the list
        clsSizes.push_back(IntInt(clsList[i]->getClusterID(),
                                  clsList[i]->getMembers().size()));
    }
    // Now sort the list based on cluster size.
    std::sort(clsSizes.begin(), clsSizes.end(),
              [](const IntInt& cls1, const IntInt& cls2) {
                  return (cls1.second < cls2.second);
              });
    // Return the sorted list of cluster sizes.
    return clsSizes;
}

#endif
