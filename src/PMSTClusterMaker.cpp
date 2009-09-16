#ifndef PMST_CLUSTER_MAKER_CPP
#define PMST_CLUSTER_MAKER_CPP

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

#include "PMSTClusterMaker.h"
#include "ESTAnalyzer.h"
#include "MSTMultiListCache.h"
#include "MSTHeapCache.h"
#include "PMSTHeapCache.h"
#include "EST.h"
#include "MPIStats.h"
#include "HeuristicChain.h"
#include "Heuristic.h"
#include "BipartiteData.h"
#include <fstream>
#include <sstream>
#include <list>
#include <math.h>

// A define to remove magic 0 (zero) in code
#define MANAGER_RANK 0

// Another define to remove magic 0 (zero) success checks.
#define NO_ERROR 0

// The default cache type to be used
char PDefCacheType[10] = "heap";

// The static instance variables for command line arguments.
int    PMSTClusterMaker::cacheSize     = 128;
bool   PMSTClusterMaker::noCacheRepop  = true;
bool   PMSTClusterMaker::strictOrder   = false;
char*  PMSTClusterMaker::inputMSTFile  = NULL;
char*  PMSTClusterMaker::outputMSTFile = NULL;
bool   PMSTClusterMaker::dontCluster   = false;
bool   PMSTClusterMaker::prettyPrint   = false;
double PMSTClusterMaker::percentile    = 1.0;
int    PMSTClusterMaker::maxUse        = 0;
char*  PMSTClusterMaker::cacheType     = PDefCacheType;

// The common set of arguments for all FW EST analyzers
arg_parser::arg_record PMSTClusterMaker::argsList[] = {
    {"--cache", "#similarity metrics to cache per EST",
     &PMSTClusterMaker::cacheSize, arg_parser::INTEGER},
    {"--no-cache-repop", "Suppress EST cache repopulation",
     &PMSTClusterMaker::noCacheRepop, arg_parser::BOOLEAN},    
    {"--percentile", "Percentile deviation to use to compute threshold value",
     &PMSTClusterMaker::percentile, arg_parser::DOUBLE},
    {"--no-order", "Disable strict order of processing messages",
     &PMSTClusterMaker::strictOrder, arg_parser::BOOLEAN},
    {"--input-mst-file", "Read MST data from file (skip parallel MST building)",
     &PMSTClusterMaker::inputMSTFile, arg_parser::STRING},
    {"--output-mst-file", "Output MST data to file",
     &PMSTClusterMaker::outputMSTFile, arg_parser::STRING},
    {"--dont-cluster", "Just generate MST data. Don't do clustering",
     &PMSTClusterMaker::dontCluster, arg_parser::BOOLEAN},
    {"--pretty-print", "Print a pretty cluster tree.",
     &PMSTClusterMaker::prettyPrint, arg_parser::BOOLEAN},
    {"--maxUse", "Set a threshold to aggressively use metrics (default=0)",
     &PMSTClusterMaker::maxUse, arg_parser::INTEGER},
    {"--cacheType", "Set type of cache (heap or mlist) to use (default=heap)",
     &PMSTClusterMaker::cacheType, arg_parser::STRING},   
    {NULL, NULL}
};

PMSTClusterMaker::PMSTClusterMaker(ESTAnalyzer *analyzer,
                                 const int refESTidx,
                                 const std::string& outputFile)
    : ClusterMaker("pmst", analyzer, refESTidx, outputFile), mst(NULL) {
    // Nothing else to be done for now.
}

PMSTClusterMaker::~PMSTClusterMaker() {
    if (mst != NULL) {
        delete mst;
    }
}

void
PMSTClusterMaker::showArguments(std::ostream& os) {
    ClusterMaker::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    os << "Options for " << name << " are:\n";
    arg_parser ap(PMSTClusterMaker::argsList);
    os << ap;
}

bool
PMSTClusterMaker::parseArguments(int& argc, char **argv) {
    arg_parser ap(PMSTClusterMaker::argsList);
    ap.check_args(argc, argv, false);
    // Ensure the cache size is at least 1.
    if (cacheSize < 1) {
        std::cerr << "Invalid cache size (must be greater than zero)\n";
        return false;
    }
    // Ensure interpretation of stringOrder flag is consistent.
    strictOrder = !strictOrder;
    // Now let the base class do processing and return the result.
    return ClusterMaker::parseArguments(argc, argv);
}

// This method is invoked only on the manager process
void
PMSTClusterMaker::estAdded(const int estIdx, std::vector<int>& repopulateList) {
    // Distribute the newly added mst node to all the workers.
    sendToWorkers(estIdx, ADD_EST);
    // A new est node has been added.  First prune our caches and
    // obtain list of nodes to be computed.
    cache->pruneCaches(estIdx, repopulateList, false);
    // Obtain and process requests to repopulate the cache from every
    // worker.
    const int MyRank = MPI_GET_RANK();
    const int WorkerCount = pData->getWorkerCount();
    for(int workerID = MyRank+1; (workerID <= MyRank+WorkerCount);
        workerID++) {
        // Wait for a message to be received from a worker and obtain
        // some status information regarding the message.
        const int sourceRank = (strictOrder ? workerID : MPI_ANY_SOURCE);
        MPI_STATUS msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        // OK, we have a valid repopulation request pending from some
        // worker. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI_INT);
        int *requestData = new int[dataSize];
        MPI_RECV(requestData, dataSize, MPI_INT,
                 msgInfo.Get_source(), REPOPULATE_REQUEST);
        // Add the population request to our repopulate vector.
        if (requestData[0] > 0) {
            std::copy(requestData + 1, requestData + dataSize - 1,
                      std::back_inserter(repopulateList));
            
        }
        // Free up memory..
        delete [] requestData;
    }
}

int
PMSTClusterMaker::managerUpdateCaches(int estIdx, const bool refreshEST) {
    // This vector is used to hold the index of ESTs whose neighbors
    // need to be recomputed as the cache for them is empty.
    std::vector<int> repopulateList;
    if (estIdx != -1) {
        estAdded(estIdx, repopulateList);
    }
    // Now that we have collected cache re-population requests from
    // all workers, add the newly added node's cache to be repopulated
    // as well, as needed.
    if (refreshEST) {
        repopulateList.push_back(estIdx);
    }
    // Now get each entry in the repopulateList to be processed.
    for(std::vector<int>::iterator curr = repopulateList.begin();
        (curr != repopulateList.end()); curr++) {
        // Send request to compute similarity list to all workers
        sendToWorkers(*curr, COMPUTE_SIMILARITY_REQUEST);
        // Participate in the cache poulation locally as well..
        populateCache(*curr);
        // Wait for the distributed cache population to complete
        // assuming the owner of this est is not the MANAGER itself.
        const int ownerRank = getOwnerProcess(*curr);
        if (ownerRank != MPI_GET_RANK()) {
            int dummy;
            TRACK_IDLE_TIME(MPI_RECV(&dummy, 1, MPI_INT, ownerRank,
                                     SIMILARITY_COMPUTATION_DONE));
        }
    }
    // Everything went on fine.
    return 0;
}

void
PMSTClusterMaker::computeNextESTidx(int& parentESTidx, int& estToAdd,
                                   float& similarity, int& alignmentData)const {
    // First send request to compute local best EST values to each
    // worker process.
    sendToWorkers(-1, COMPUTE_MAX_SIMILARITY_REQUEST);
    // Now compute local best similarity
    cache->getBestEntry(parentESTidx, estToAdd, similarity, alignmentData);
    // Receive similarity entry from
    const int MyRank = MPI_GET_RANK();
    const int WorkerCount = pData->getWorkerCount();
    
    for(int rank = MyRank+1; (rank <= MyRank+WorkerCount); rank++) {
        // Choose worker rank depending on strict ordering scheme..
        const int workerRank = (strictOrder ? rank : MPI_ANY_SOURCE);
        // Get the local simlarity information from another worker.
        int remoteData[4];
        TRACK_IDLE_TIME(MPI_RECV(remoteData, 4, MPI_INT,
                                 workerRank, MAX_SIMILARITY_RESPONSE));
        // Undo the fudge on similarity done at the sender end.
        const float remoteSim = *((float *) (remoteData + 2));
        // Use a better or first valid entry
        if (analyzer->compareMetrics(remoteSim, similarity) ||
            (estToAdd == -1)) {
            // Found a higher similarity or a shorter distance in a
            // remote process or this is the first valid entry thusfar
            similarity   = remoteSim;
            parentESTidx = remoteData[0];
            estToAdd     = remoteData[1];
            alignmentData= remoteData[3];
        }
    } 
}

void
PMSTClusterMaker::addMoreChildESTs(const int parentESTidx, int& estToAdd,
                                  float &metric, int& alignmentData,
                                  int& pendingESTs) {
    // variable to track parent of next EST to be added.
    int newParent = -1;
    
    do {
        // Add the current node to the MST.
        ASSERT ( estToAdd != -1 );
        mst->addNode(parentESTidx, estToAdd, metric, alignmentData);
        if (--pendingESTs == 0) {
            // All EST's have been added to the MST.  Nothing more to
            // do.  So break out of the while loop.
            break;
        }
        // Have the manager and worker update their caches to
        // appropriately prune entries for the newly added EST.
        managerUpdateCaches(estToAdd, false);
        // Now determine if there is another child that can be added using
        // an helper method and without destroying current values.
        int newChild       = -1;
        float childMetric  = 0;
        int childAlignment = 0;
        computeNextESTidx(newParent, newChild, childMetric, childAlignment);
        // Check if we have a child with a useful alignment data.
        if ((newParent == parentESTidx) &&
            (analyzer->compareMetrics(childMetric, (float) maxUse))) {
            // Found a child of the same parent with a good metric. So
            // update parameters with the new child.
            estToAdd      = newChild;
            metric        = childMetric;
            alignmentData = childAlignment;
        } else {
            // No more children to add!
            newParent = -1;
        }
    } while (newParent == parentESTidx);
}

int
PMSTClusterMaker::manager() {
    // The number of pending nodes to be added to the MST.
    int pendingESTs = pData->estCount;
    // The minimum spanning tree that is built by this manager.
    int dummy;
    mst = new MST(pendingESTs, analyzer->getAlignmentData(dummy));
    // For building sub-MSTs, we make the starting EST the reference EST.
    int parentESTidx    = -1;
    int estToAdd        = pData->startESTidx;
    float metric        = 0;
    int   alignmentInfo = 0;
    do {
        // Add the EST to the MST vector, if needed.
        if ((maxUse == -1) || (parentESTidx == -1)) {
            ASSERT ( estToAdd != -1 );
            mst->addNode(parentESTidx, estToAdd, metric, alignmentInfo);
            if (--pendingESTs == 0) {
                // All EST's have been added to the MST.  Nothing more to
                // do.  So break out of the while loop.
                break;
            }
        }
        // process requests to repopulate caches and populate cache
        // for the newly addded EST entry.
        managerUpdateCaches(estToAdd);
        // Now broadcast request to all workers to provide their best
        // choice for the next EST id to be added to the MST using a
        // helper method.
        computeNextESTidx(parentESTidx, estToAdd, metric, alignmentInfo);
        ASSERT( parentESTidx != -1 );
        ASSERT( estToAdd     != -1 );
        if (maxUse != -1) {
            // Try to add as many ESTs as possible rooted at the given
            // parentESTidx using a helper method.
            addMoreChildESTs(parentESTidx, estToAdd, metric,
                             alignmentInfo, pendingESTs);
        }
    } while (pendingESTs > 0);
    
    // Broad cast an estIdx of -1 to all the workers to indicate that
    // MST building is done.
    sendToWorkers(-1, ADD_EST);
    printf("Complete %d\n", MPI_GET_RANK());
    // All done with no problems...
    return NO_ERROR;
}

int
PMSTClusterMaker::worker() {
    // In the first step the worker first receives the index of the
    // EST was just added to the MST.  If the estNode added is -1,
    // then the MST building process is done.
    int estAdded = -1;
    // Wait for the Manager to send requests to this worker to perform
    // different tasks.
    MPI_STATUS msgInfo;
    const int LocalManagerRank = pData->getPartitionManager();
    do {
        // Wait for manager to send us a work request.  Since we are
        // waiting it should be tracked under idle time.
        MPI_PROBE(LocalManagerRank, MPI_ANY_TAG, msgInfo);
        if (msgInfo.Get_tag() == COMPUTE_SIMILARITY_REQUEST) {
            // Read the actual message first.  Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            int estIdx = -1;
            MPI_RECV(&estIdx, 1, MPI_INT, LocalManagerRank,
                     COMPUTE_SIMILARITY_REQUEST);
            // Perform the necessary operations.
            populateCache(estIdx);
        } else if (msgInfo.Get_tag() == COMPUTE_MAX_SIMILARITY_REQUEST) {
            // Read the actual message first. Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            int dummy = 0;
            MPI_RECV(&dummy, 1, MPI_INT, LocalManagerRank,
                     COMPUTE_MAX_SIMILARITY_REQUEST);
            int   bestEntry[4];
            float similarity = 0;
            // Get the best possible local similarity match.
            cache->getBestEntry(bestEntry[0], bestEntry[1],
                                similarity, bestEntry[3]);
            // Fudge the similarity into the bestEntry array to
            // transmitt the necessary information to the manager.
            // Maybe there is a cleaner way to do it too...
            int *temp    = reinterpret_cast<int*>(&similarity);
            bestEntry[2] = *temp;
            MPI_SEND(bestEntry, 4, MPI_INT, LocalManagerRank,
                     MAX_SIMILARITY_RESPONSE);
        } else if (msgInfo.Get_tag() == ADD_EST) {
            // The manager has broad casted the next est to be added.
            MPI_RECV(&estAdded, 1, MPI_INT, LocalManagerRank, ADD_EST);
            if (estAdded == -1) {
                // No more ESTs to add.  Clustering is done.  So it is time
                // for this worker to stop too.
                break;
            }
            // A new est node has been added.  First prune our caches and post
            // any new cache computation requests back to the Manager.
            std::vector<int> repopulateList;
            cache->pruneCaches(estAdded, repopulateList);
            // Send the repopulate list to the manager.
            MPI_SEND(&repopulateList[0], repopulateList.size(),
                     MPI_INT, LocalManagerRank, REPOPULATE_REQUEST);
        }
    } while (estAdded != -1);
    // Everything went on without a hitch.
    return NO_ERROR;
}

int
PMSTClusterMaker::getOwnerProcess(const int estIdx) const {
    const int MyRank = MPI_GET_RANK();
    const int ProcessCount = pData->getWorkerCount() + 1;
    if (ProcessCount > 1) {
        // This is a bipartite graph with multiple processes
        // ESTs in each partition are split up equally among the processes
        BipartiteData* bpData = static_cast<BipartiteData*> (pData);
        int startIdx1 = bpData->partition1->startESTidx;
        int startIdx2 = bpData->partition2->startESTidx;
        int rank = 0;
        if (estIdx < startIdx2) {
            int factor = bpData->partition1->estCount / ProcessCount;
            int rem = bpData->partition1->estCount % ProcessCount;
            rank = (estIdx - startIdx1 - rem) / factor;
            if (rank < 0) rank = 0;
            //rank = (estIdx - startIdx1) / ((bpData->partition1->estCount
            //                             / ProcessCount)+1);
        } else {
            int factor = bpData->partition2->estCount / ProcessCount;
            int rem = bpData->partition2->estCount % ProcessCount;
            rank = (estIdx - startIdx2 - rem) / factor;
            if (rank < 0) rank = 0;
            //rank = (estIdx - startIdx2) / ((bpData->partition2->estCount
            //                             / ProcessCount)+1);
        }
        return rank + bpData->getPartitionManager();
    }
    return MyRank;
}

float
PMSTClusterMaker::analyze(const int otherEST) {
    return analyzer->analyze(otherEST);
}

void
PMSTClusterMaker::getOwnedESTidx(const int estIdx, int& startIndex,
                                 int& endIndex) {
    if (pData->getPartitionManager() == -1) {
        // Normal partition, this process owns all ESTs
        startIndex = pData->startESTidx;
        endIndex = pData->startESTidx + pData->estCount;
    } else {
        // This is a bipartite graph with multiple processes
        // ESTs in each partition are split up equally among the processes
        // Find process count and relative rank
        const int ProcessCount = pData->getWorkerCount() + 1;
        int rank = MPI_GET_RANK() - pData->getPartitionManager();
        BipartiteData* bpData = static_cast<BipartiteData*> (pData);
        int startIdx1 = bpData->partition1->startESTidx;
        int startIdx2 = bpData->partition2->startESTidx;
        if (estIdx < startIdx2) {
            // EST is in first partition, so we analyze the second partition
            int factor = bpData->partition2->estCount / ProcessCount;
            int rem = bpData->partition2->estCount % ProcessCount;
            startIndex = startIdx2 + (factor * rank);
            endIndex = startIndex + factor;
            if (rank == 0) endIndex += rem;
            else {
                startIndex += rem;
                endIndex += rem;
            }
        } else {
            // EST is in second partition, so we analyze the first partition
            int factor = bpData->partition1->estCount / ProcessCount;
            int rem = bpData->partition1->estCount % ProcessCount;
            startIndex = startIdx1 + (factor * rank);
            endIndex = startIndex + factor;
            if (rank == 0) endIndex += rem;
            else {
                startIndex += rem;
                endIndex += rem;
            }
        }
        //printf("%d %d %d\n", estIdx, startIndex, endIndex);
    }
}

void
PMSTClusterMaker::populateCache(const int estIdx, SMList* metricList) {
    // First determine the list of ESTs that this process must deal
    // with using the helper method.
    int startESTidx, endESTidx;
    getOwnedESTidx(estIdx, startESTidx, endESTidx);
    // Setup the reference estIdx in the analyzer which given the
    // analyzer a chance to optimize initialization.
    analyzer->setReferenceEST(estIdx);
    // Now compute similarity metric and store information in a SMList
    // data structure.
    SMList smList;
    const  float InvalidMetric = analyzer->getInvalidMetric();
    // Flag to ensure only one invalid metric gets added
    bool needInvalidMetric = true;
    
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if ((otherIdx == estIdx) || (cache->isESTinMST(otherIdx))) {
            // This EST entry can be ignored as this similarity metric
            // is not needed for MST construction.
            continue;
        }
        // Get similarity/distance metric.
        const float metric = analyze(otherIdx);
        ASSERT ( metric >= 0 );
        // Obtain any alignemnt data that the analyzer may have.
        int alignmentData = 0;
        analyzer->getAlignmentData(alignmentData);
        // Add only the first invalid entry. One is enough to do build
        // a valid MST. There is no need for multiple vestigial
        // entries.
        const bool isMetricOK= analyzer->compareMetrics(metric, InvalidMetric);
        if ((isMetricOK) || (needInvalidMetric)) {
            // Add the information to metric list
            smList.push_back(CachedESTInfo(estIdx, otherIdx,
                                           metric, alignmentData));
            // If this is an invalid entry, then one is enough. Don't
            // add more.
            needInvalidMetric = (pData->getPartitionManager() != -1) ||
                (needInvalidMetric && isMetricOK);                
        }
    }
    // Preprocess the SMList to make it optimal for the MSTCache to
    // process in a distributed manner.
    cache->preprocess(smList);
    // Now further process smList...
    const int ownerRank = getOwnerProcess(estIdx);
    if (ownerRank != MPI_GET_RANK()) {
        // This process is not the owner.  In this case, just send the
        // smList to the remote owner process.  There is some fudging
        // of data types going on here using the assumption that
        // sizeof(int) == sizeof(float).  Maybe can make this code
        // more MPI friendly.
        
        // First ensure there is at least one entry to work with.
        if (smList.empty()) {
            smList.push_back(CachedESTInfo(-1, -1, -1.0f, -1));
        }
        MPI_SEND(&smList[0], smList.size() * sizeof(CachedESTInfo),
                 MPI_CHAR, ownerRank, SIMILARITY_LIST);
        // Nothing futher to do if the process is not the owner for
        // this EST.
        return;
    }
    // When control drops here we are dealing with the process that
    // owns this EST.  This process must obtain similarity_lists from
    // other processes and merge them together into its local cache.
    // First merge the locally computed smList first.
    cache->mergeList(estIdx, smList);
    // Add the entries to the metricList if requested
    if (metricList != NULL) {
        metricList->insert(metricList->end(), smList.begin(), smList.end());
    }
    // Add data to the parameter
    // Obtain similarity lists from all other processes (other than
    // ourselves).

    const int LocalManagerRank = pData->getPartitionManager();
    if (LocalManagerRank == -1) {
        // Serial computation of this MST, so we're finished here
        return;
    }
    const int MyRank = MPI_GET_RANK();
    for(int pid = LocalManagerRank;
        (pid <= LocalManagerRank+pData->getWorkerCount()); pid++) {
        if (pid == MyRank) {
            // Can't get a message from ourselves...
            continue;
        }
        // Choose actual rank depending on strict ordering scheme..
        const int rank = (strictOrder ? pid : MPI_ANY_SOURCE);
        MPI_STATUS msgInfo;
        // Wait to receive similarity list.  Since we are waiting it
        // should be logged as idle time for this process.
        MPI_PROBE(rank, SIMILARITY_LIST, msgInfo);
        // OK, we have a valid similarity list pending from some other
        // process. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI_CHAR);
        SMList remoteList(dataSize / sizeof(CachedESTInfo));
        // The following call is a kludge with MPI/STL data types
        // based on several language assumptions.  This part could be
        // cleaned up to be more portable later on.
        MPI_RECV(&remoteList[0], dataSize, MPI_CHAR,
                 msgInfo.Get_source(), SIMILARITY_LIST);
        // Merge the list we got from the remote process with our
        // local cache information if the list has a valid entry.
        if (hasValidSMEntry(remoteList)) {
            cache->mergeList(estIdx, remoteList);
            // Add the entries to the metricList if requested
            if (metricList != NULL) {
                metricList->insert(metricList->end(), remoteList.begin(),
                                   remoteList.end());
            }
        }
    }
    // Now finally let the manager know that the round of similarity
    // computation is all completed (assuming this process itself is
    // not the manager).
    if (MPI_GET_RANK() != LocalManagerRank) {
        const int dummy = -1;
        MPI_SEND(&dummy, 1, MPI_INT, LocalManagerRank,
                 SIMILARITY_COMPUTATION_DONE);
    }
}

void
PMSTClusterMaker::getOwnedPartition() {
    const int CommSize         = MPI_GET_SIZE();
    // Determine the number of partitions.
    // The EST set is partitioned into a number of roughly equal subsets
    // based on the number of processes we have available.
    // Based on the number of partitions s, we will also have a total of
    // s(s-1)/2 bipartite graphs.  Since the number of analyses performed
    // in MST construction on the bipartite graphs is greater, we would like
    // to assign 2 processes to work on each bipartite graph if possible.
    // Thus the formula: (this formula will be optimized after testing)
    const int NumPartitions    = sqrt(CommSize); // rounds down
    const int MyRank           = MPI_GET_RANK();    
    const int ESTsPerPartition = EST::getESTList().size() / NumPartitions;
    const int ExtraESTs        = EST::getESTList().size() % NumPartitions;
    int startESTidx = 0;
    int estCount = 0;

    if (MyRank < NumPartitions) {
        // Assign to a normal partition
        if (MyRank <= ExtraESTs) {
            // The previous processes have one extra ESTs as number of
            // ESTs are not evenly divisible by number of processes.  So
            // account for this.
            startESTidx = ((ESTsPerPartition + 1) * MyRank);
        } else {
            startESTidx = MyRank * ESTsPerPartition + ExtraESTs;
        }
        estCount = ESTsPerPartition;
        if (MyRank < ExtraESTs) {
            // This process owns one extra EST as ESTs are not evenly
            // divisible by the number of processes.
            estCount++;
        }
        pData = new PartitionData(startESTidx, estCount);
    } else {
        // Handle bipartite assignment.  Need to let managers know how
        // many workers they will have and make sure workers know their
        // manager.  Manager can figure out worker ranks from the
        // worker count and their rank.
        int bpCount = NumPartitions * (NumPartitions-1) / 2;
        int bpProcessCount = CommSize - NumPartitions;
        int procsPerBP = bpProcessCount / bpCount;
        int extraProcs = bpProcessCount % bpCount;
        int rank = MyRank - NumPartitions;

        // First account for extra processes
        if (rank < extraProcs*(procsPerBP+1)) {
            procsPerBP += 1;
        }
        int manager = rank;
        int workerCount = procsPerBP - 1;
        if ((rank % procsPerBP) != 0) {
            // Worker for bipartite
            manager = rank - (rank % procsPerBP);
        }
        // Figure out which graphs make up the bipartite
        int bpID = 0;
        int id1 = 0;
        int id2 = 0;
        if (rank >= extraProcs) {
            bpID = (rank - extraProcs) / procsPerBP; // not sure about this
        } else {
            bpID = rank / procsPerBP;
        }
        id1 = ((bpID) * 2) / NumPartitions;
        if (bpID < NumPartitions-1) {
            id2 = bpID+1;
        } else {
            id2 = ((bpID-1) % (NumPartitions-id1)) + id1;
        }
           
        //printf("%d %d %d %d\n", NumPartitions, bpID, id1, id2);
        
        // Assign to a normal partition
        if (id1 <= ExtraESTs) {
            // The previous processes have one extra ESTs as number of
            // ESTs are not evenly divisible by number of processes.  So
            // account for this.
            startESTidx = ((ESTsPerPartition + 1) * id1);
        } else {
            startESTidx = id1 * ESTsPerPartition + ExtraESTs;
        }
        estCount = ESTsPerPartition;
        if (id1 < ExtraESTs) {
            // This process owns one extra EST as ESTs are not evenly
            // divisible by the number of processes.
            estCount++;
        }

        int startESTidx2 = 0;
        int estCount2 = 0;
        if (id2 <= ExtraESTs) {
            // The previous processes have one extra ESTs as number of
            // ESTs are not evenly divisible by number of processes.  So
            // account for this.
            startESTidx2 = ((ESTsPerPartition + 1) * id2);
        } else {
            startESTidx2 = id2 * ESTsPerPartition + ExtraESTs;
        }
        estCount2 = ESTsPerPartition;
        if (id2 < ExtraESTs) {
            // This process owns one extra EST as ESTs are not evenly
            // divisible by the number of processes.
            estCount2++;
        }

        // Assign to the bipartite graph
        pData = new BipartiteData(startESTidx, estCount+estCount2,
                                  new PartitionData(startESTidx, estCount),
                                  new PartitionData(startESTidx2, estCount2),
                                  manager+NumPartitions, workerCount);
        
        //printf("%d %d %d %d\n", startESTidx, estCount, startESTidx2, estCount2);
    }
}

void
PMSTClusterMaker::displayStats(std::ostream& os) {
    // Dump cache usage statistics for this process.
    cache->displayStats(os, MPI_GET_RANK());
    // Display MPI usage statistics.
    MPIStats::displayStats(os);
}

int
PMSTClusterMaker::mergeManager(MSTCluster& rootCluster, const int threshold) {
    // Create a new cache for constructing the "full" MST.
    // This cache will simply hold all edges from the merged MSTs.
    PMSTHeapCache* pCache
        = new PMSTHeapCache(EST::getESTList().size(), 0,
                            EST::getESTList().size(), analyzer);
    // Extract data from our MST and put into a new SMList.

    NodeList nodes = mst->getNodes();
    SMList smList;

    // Iterate through the nodelist
    for(NodeList::const_iterator entry = nodes.begin();
        (entry != nodes.end()); entry++) {
        if (entry->parentIdx != -1) {
            // Add the information to metric list
            smList.push_back(CachedESTInfo(entry->parentIdx, entry->estIdx,
                             entry->metric, entry->alignmentMetric));
        }
    }    

    // Old MST is no longer needed.
    delete mst;

    // Merge the SMList with the cache.
    
    pCache->mergeList(0, smList);

    // Receive all SMLists from workers and merge them with cache.
    // Need to calculate number of subgraphs first.
    int numPartitions = sqrt(MPI_GET_SIZE());
    int numGraphs = numPartitions + ((numPartitions * (numPartitions-1)) / 2);
    for (int count = 1; count < numGraphs; count++) {
        MPI_STATUS msgInfo;
        // Wait to receive similarity list.  Since we are waiting it
        // should be logged as idle time for this process.
        MPI_PROBE(MPI_ANY_SOURCE, SIMILARITY_LIST, msgInfo);
        // OK, we have a valid similarity list pending from some other
        // process. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI_CHAR);
        SMList remoteList(dataSize / sizeof(CachedESTInfo));
        // The following call is a kludge with MPI/STL data types
        // based on several language assumptions.  This part could be
        // cleaned up to be more portable later on.
        MPI_RECV(&remoteList[0], dataSize, MPI_CHAR,
                 msgInfo.Get_source(), SIMILARITY_LIST);
        // Merge the list we got from the remote process with our
        // local cache information.
        pCache->mergeList(0, remoteList);
    }

    // Now construct the MST using the new cache.

    // The number of pending nodes to be added to the MST.
    int pendingESTs = EST::getESTList().size();

    // union-find stuff
    int* parent = new int[pendingESTs];
    bool* root  = new bool[pendingESTs];
    for (int i = 0; i < pendingESTs; i++) {
        root[i] = true;
        parent[i] = 1;
    }

    // The minimum spanning tree that is built by this manager.
    int dummy;
    mst = new MST(pendingESTs, analyzer->getAlignmentData(dummy));

    int parentESTidx    = -1;
    int estToAdd        = refESTidx;
    float metric        = 0;
    int   alignmentInfo = 0;

    pCache->popBestEntry(parentESTidx, estToAdd, metric, alignmentInfo);
    //pCache->pruneMergedCaches(parentESTidx, estToAdd);
    // Normally would check to prevent cycles, but we're guaranteed this
    // is a valid entry since it is the first entry
    
    // Add the first parent as the root of the MST with metric of 0
    mst->addNode(-1, parentESTidx, 0, 0);
    pendingESTs--;

    int i = parentESTidx;
    int j = estToAdd;

    bool reachedThreshold = false;

    do {
        if (!reachedThreshold && metric > threshold) {
            // Hit the threshold.  All edges from here out will be greater
            // (as guaranteed by Kruskal's algorithm)
            reachedThreshold = true;
            // So, we should make the clusters now
            rootCluster.makeMergedClusters(EST::getESTList().size(),
                                           parent, root);
        }
        // Add the edge
        mst->addNode(parentESTidx, estToAdd, metric, alignmentInfo);
        if (--pendingESTs == 0) {
            // All EST's have been added to the MST.  Nothing more to
            // do.  So break out of the while loop.
            break;
        }

        // Perform the union operation
        if (parent[i] < parent[j]) { // compare sizes
            parent[j] += parent[i];
            root[i] = false;
            parent[i] = j;
        } else {
            parent[i] += parent[j];
            root[j] = false;
            parent[j] = i;
        }
        
        //pCache->pruneMergedCaches(parentESTidx, estToAdd);

        do {
            // Compute local best similarity, for next EST id to be added
            pCache->popBestEntry(parentESTidx, estToAdd, metric,
                                     alignmentInfo);
            i = parentESTidx;
            j = estToAdd;
            int itemp = parentESTidx;
            int jtemp = estToAdd;
            while (!root[i]) {
                i = parent[i];
            }
            // path compression
            while (itemp != i) {
                int next = parent[itemp];
                parent[itemp] = i;
                itemp = next;
            }
            while (!root[j]) {
                j = parent[j];
            }
            // path compression
            while (jtemp != j) {
                int next = parent[jtemp];
                parent[jtemp] = j;
                jtemp = next;
            }
        } while (i == j);
      
        ASSERT( parentESTidx != -1 );
        ASSERT( estToAdd     != -1 );
        
    } while (pendingESTs > 0);
    
    // Clean up.
    delete [] parent;
    delete [] root;
    delete pCache;
    
    // All done with no problems...
    return NO_ERROR;
}

int
PMSTClusterMaker::mergeWorker() {
    // Extract data from our MST and put into a new SMList.
    NodeList nodes = mst->getNodes();
    SMList smList;

    // Iterate through the nodelist
    for(NodeList::const_iterator entry = nodes.begin();
        (entry != nodes.end()); entry++) {
        if (entry->parentIdx != -1) {
            // Add the information to metric list
            smList.push_back(CachedESTInfo(entry->parentIdx, entry->estIdx,
                             entry->metric, entry->alignmentMetric));
        }
    }

    // Send the new SMList to the manager.
    MPI_SEND(&smList[0], smList.size() * sizeof(CachedESTInfo),
             MPI_CHAR, MANAGER_RANK, SIMILARITY_LIST);

    // Delete the MST.
    delete mst;
    mst = NULL;

    // Everything went on without a hitch.
    return NO_ERROR;
}

int
PMSTClusterMaker::makeClusters() {
    int result = NO_ERROR;

    MSTCluster root; // make this parameterizable

    // All the parallel MPI processes load the data for processing
    // into memory.
    if ((result = analyzer->initialize()) != NO_ERROR) {
        // Error occured during initialization. Bail out.
        return result;
    }
    
    int totalSuccesses = 0;

    if (inputMSTFile == NULL) {
        // No input MST file supplied.  That means the MST must be
        // built. The local cache that contains information to build
        // the MST.  The cache for that contains information regarding
        // ESTs under the ownership of this process.  The cache
        // facilitates rapid construction of ESTs but minimizing the
        // number of times similarity metrics need to computed.  The
        // cache needs to know the starting indedx of the EST whose
        // data must be cached by this process.  So this information
        // is computed first.

        // Processes get the ESTs for which they are responsible.
        getOwnedPartition();
        if (!strcmp(cacheType, "mlist")) {
            cache = new MSTMultiListCache(EST::getESTList().size(),
                                          pData->startESTidx,
                                          pData->ownedESTCount, analyzer,
                                          !noCacheRepop, cacheSize);
        } else {
            cache = new MSTHeapCache(EST::getESTList().size(),
                                     pData->startESTidx,
                                     pData->ownedESTCount, analyzer,
                                     !noCacheRepop, cacheSize);
        }
        
        // Processes construct MSTs on their subgraphs in parallel.

        if (pData->getPartitionManager() == -1 ||
            pData->getPartitionManager() == MPI_GET_RANK()) {
            //printf("mgr %d %d\n", pData->startESTidx, pData->estCount);
            result = manager();
        } else {
            //printf("wkr %d %d\n", pData->startESTidx, pData->estCount);
            result = worker();
        }

        // Collaboratively calculate the clustering threshold.

        int threshold = 40; // default
        Heuristic* tv = (HeuristicChain::getHeuristicChain())
            ->getHeuristic("tv");
        if (tv != NULL) {
            // We are good to go for summing the statistics
            int tvSuccesses = tv->getSuccessCount();
            // Add the manager's successes first
            totalSuccesses = tvSuccesses;

            // Get each worker's number of successes and add them
            if (MPI_GET_RANK() == MANAGER_RANK) {
                for (int i = 1; i < MPI_GET_SIZE(); i++) {
                    TRACK_IDLE_TIME(MPI_RECV(&tvSuccesses, 1, MPI_INT,
                                             MPI_ANY_SOURCE,
                                             COMPUTE_TOTAL_ANALYSIS_COUNT));
                    totalSuccesses+=tvSuccesses;
                }
                // Calculate threshold
                threshold = root.calculateThreshold(EST::getESTList().size(),
                                        percentile, totalSuccesses, analyzer);
            } else {
                // Workers send
                MPI_SEND(&tvSuccesses, 1, MPI_INT, MANAGER_RANK,
                         COMPUTE_TOTAL_ANALYSIS_COUNT);
            }
        }

        // Display statistics by writing all the data to a string
        // stream and then flushing the stream.  This tries to working
        // around (it is not perfect solution) interspersed data from
        // multiple MPI processes

        std::ostringstream buffer;
        displayStats(buffer);
        std::cout << buffer.str() << std::endl;
        // Delete the cache as we no longer needed it.
        delete cache;

        // Processes must now send their results to the manager, so that
        // the manager can merge the subgraph MSTs and make the full MST.
        // Act as manager or worker depending on MPI rank.
        if (MPI_GET_RANK() == MANAGER_RANK) {
            // Get this MPI process to act as the manager.
            result = mergeManager(root, threshold);
        } else if (pData->getPartitionManager() == -1 ||
                   pData->getPartitionManager() == MPI_GET_RANK()) {
            // Managers of partitions must contribute their MSTs to the
            // merged MST, and act as "workers" in this process.
            result = mergeWorker();
        }
        // partition data no longer needed
        delete pData;
    } else {
        // In this case the MST is to be read from an input file.
        // So let's do that.
        mst = MST::deSerialize(inputMSTFile);
    }

    if ((result != NO_ERROR) || (mst == NULL)) {
        // Some error occured.  Can't proceed further.
        return result;
    }
    ASSERT ( mst != NULL );
    // Dump MST out to the specified output file.
    if ((outputMSTFile != NULL) && (mst != NULL)) {
        mst->serialize(outputMSTFile, (inputMSTFile != NULL) ? inputMSTFile :
                       analyzer->getInputFileName());
    }
    // Output the clusters if so desired.
    if (!dontCluster) {
        //root.makeClusters(mst->getNodes(), percentile, totalSuccesses);
                    
        // Redirect cluster output to outputFile as needed.
        std::ofstream outFile;
        if (!outputFileName.empty()) {
            outFile.open(outputFileName.c_str());
            if (!outFile.good()) {
                std::cerr << "Error opening output file " << outputFileName
                          << " for writing cluster data.\n"
                          << "Cluster data will be dumped on stdout.\n";
            }
        }
        if (prettyPrint) {
            root.printClusterTree((outFile.is_open() && outFile.good()) ?
                                  outFile : std::cout);
        } else {
            // No pretty printing. Just dump the info out.
            ((outFile.is_open() && outFile.good()) ? outFile : std::cout)
                << root << std::endl;
        }
    }
    // Return result from behaving as a manager or a worker.
    return result;
}

void
PMSTClusterMaker::sendToWorkers(int data, const int tag) const {
    if (pData->getPartitionManager() == MPI_GET_RANK()) {
        const int MyRank = MPI_GET_RANK();
        for (int i = MyRank+1; i <= MyRank+pData->getWorkerCount(); i++) {
            MPI_SEND(&data, 1, MPI_INT, i, tag);
        }
    }
}

bool
PMSTClusterMaker::hasValidSMEntry(const SMList& list) const {
    if (list.empty()) {
        // If the list is empty, then there is nothing further to be
        // done.
        return false;
    }
    if (list.size() > 1) {
        // If the list has more than one entry, then the list must
        // have at least one valid entry.
        return true;
    }
    if (list[0].estIdx == -1) {
        // The list has 1 entry and the first entry's index is -1 that
        // means this list does not have a even one valid entry.
        return false;
    }
    // The list does have one valid entry.
    return true;
}

#endif
