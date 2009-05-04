#ifndef MST_CLUSTER_MAKER_CPP
#define MST_CLUSTER_MAKER_CPP

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

#include "MSTClusterMaker.h"
#include "ESTAnalyzer.h"
#include "MSTCluster.h"
#include "EST.h"
#include "MPIStats.h"
#include <fstream>

// A define to remove magic 0 (zero) in code
#define MANAGER_RANK 0

// Another define to remove magic 0 (zero) success checks.
#define NO_ERROR 0

// The static instance variables for command line arguments.
int    MSTClusterMaker::cacheSize     = 128;
bool   MSTClusterMaker::noCacheRepop  = true;
bool   MSTClusterMaker::strictOrder   = false;
char*  MSTClusterMaker::inputMSTFile  = NULL;
char*  MSTClusterMaker::outputMSTFile = NULL;
bool   MSTClusterMaker::dontCluster   = false;
bool   MSTClusterMaker::prettyPrint   = false;
double MSTClusterMaker::percentile    = 1.0;
int    MSTClusterMaker::maxUse        = 1;

// The common set of arguments for all FW EST analyzers
arg_parser::arg_record MSTClusterMaker::argsList[] = {
    {"--cache", "#similarity metrics to cache per EST",
     &MSTClusterMaker::cacheSize, arg_parser::INTEGER},
    {"--no-cache-repop", "Suppress EST cache repopulation",
     &MSTClusterMaker::noCacheRepop, arg_parser::BOOLEAN},    
    {"--percentile", "Percentile deviation to use to compute threshold value",
     &MSTClusterMaker::percentile, arg_parser::DOUBLE},
    {"--no-order", "Disable strict order of processing messages",
     &MSTClusterMaker::strictOrder, arg_parser::BOOLEAN},
    {"--input-mst-file", "Read MST data from file (skip parallel MST building)",
     &MSTClusterMaker::inputMSTFile, arg_parser::STRING},
    {"--output-mst-file", "Output MST data to file",
     &MSTClusterMaker::outputMSTFile, arg_parser::STRING},
    {"--dont-cluster", "Just generate MST data. Don't do clustering",
     &MSTClusterMaker::dontCluster, arg_parser::BOOLEAN},
    {"--pretty-print", "Print a pretty cluster tree.",
     &MSTClusterMaker::prettyPrint, arg_parser::BOOLEAN},
    {"--maxUse", "Set a threshold to aggressively use metrics",
     &MSTClusterMaker::maxUse, arg_parser::INTEGER},    
    {NULL, NULL}
};

MSTClusterMaker::MSTClusterMaker(ESTAnalyzer *analyzer,
                                 const int refESTidx,
                                 const std::string& outputFile)
    : ClusterMaker("mst", analyzer, refESTidx, outputFile), mst(NULL) {
    // Nothing else to be done for now.
}

MSTClusterMaker::~MSTClusterMaker() {
    if (mst != NULL) {
        delete mst;
    }
}

void
MSTClusterMaker::showArguments(std::ostream& os) {
    ClusterMaker::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    os << "Options for " << name << " are:\n";
    arg_parser ap(MSTClusterMaker::argsList);
    os << ap;
}

bool
MSTClusterMaker::parseArguments(int& argc, char **argv) {
    arg_parser ap(MSTClusterMaker::argsList);
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
MSTClusterMaker::estAdded(const int estIdx, std::vector<int>& repopulateList) {
    // Distribute the newly added mst node to all the workers.
    sendToWorkers(estIdx, ADD_EST);
    // A new est node has been added.  First prune our caches and
    // obtain list of nodes to be computed.
    cache->pruneCaches(estIdx, repopulateList, false);
    // Obtain and process requests to repopulate the cache from every
    // worker.
    const int WorkerCount = MPI::COMM_WORLD.Get_size();
    for(int workerID = 1; (workerID < WorkerCount); workerID++) {
        // Wait for a message to be received from a worker and obtain
        // some status information regarding the message.
        const int sourceRank = (strictOrder ? workerID : MPI_ANY_SOURCE);
        MPI::Status msgInfo;
        MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
        // OK, we have a valid repopulation request pending from some
        // worker. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI::INT);
        int *requestData = new int[dataSize];
        MPI_RECV(requestData, dataSize, MPI::INT,
                 msgInfo.Get_source(), REPOPULATE_REQUEST);
        // Add the poulation request to our repopulate vector.
        if (requestData[0] > 0) {
            std::copy(requestData + 1, requestData + dataSize - 1,
                      std::back_inserter(repopulateList));
            
        }
        // Free up memory..
        delete [] requestData;
    }
}

int
MSTClusterMaker::managerUpdateCaches(int estIdx, const bool refreshEST) {
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
        if (ownerRank != MANAGER_RANK) {
            int dummy;
            TRACK_IDLE_TIME(MPI_RECV(&dummy, 1, MPI_INT, ownerRank,
                                     SIMILARITY_COMPUTATION_DONE));
        }
    }
    // Everything went on fine.
    return MPI::SUCCESS;
}

void
MSTClusterMaker::computeNextESTidx(int& parentESTidx, int& estToAdd,
                                   float& similarity, int& alignmentData)const {
    // First send request to compute local best EST values to each
    // worker process.
    sendToWorkers(-1, COMPUTE_MAX_SIMILARITY_REQUEST);
    // Now compute local best similarity
    cache->getBestEntry(parentESTidx, estToAdd, similarity, alignmentData);
    // Receive similarity entry from 
    const int ProcessCount = MPI::COMM_WORLD.Get_size();
    for(int rank = 1; (rank < ProcessCount); rank++) {
        // Choose worker rank depending on strict ordering scheme..
        const int workerRank = (strictOrder ? rank : MPI_ANY_SOURCE);
        // Get the local simlarity information from another worker.
        int remoteData[4];
        TRACK_IDLE_TIME(MPI_RECV(remoteData, 4, MPI::INT,
                                 workerRank, MAX_SIMILARITY_RESPONSE));
        // Undo the fudge on similarity done at the sender end.
        const float remoteSim = *((float *) (remoteData + 2));
        if (analyzer->compareMetrics(remoteSim, similarity)) {
            // Found a higher similarity or a shorter distance in a
            // remote process!
            similarity   = remoteSim;
            parentESTidx = remoteData[0];
            estToAdd     = remoteData[1];
            alignmentData= remoteData[3];
        }
    }
}

void
MSTClusterMaker::addMoreChildESTs(const int parentESTidx, int& estToAdd,
                                  float &metric, int& alignmentData,
                                  int& pendingESTs) {
    // variable to track parent of next EST to be added.
    int newParent = -1;
    
    do {
        // Add the current node to the MST.
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
            (analyzer->compareMetrics(childMetric, maxUse))) {
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
MSTClusterMaker::manager() {
    // The number of pending nodes to be added to the MST.
    int pendingESTs = EST::getESTList().size();
    // The minimum spanning tree that is built by this manager.
    int dummy;
    mst = new MST(pendingESTs, analyzer->getAlignmentData(dummy));
    // Kick off all activities by adding the reference EST as the root
    // of the MST with a similarity metric of 0.
    int parentESTidx    = -1;
    int estToAdd        = refESTidx;
    float metric        = 0;
    int   alignmentInfo = 0;
    do {
        // Add the EST to the MST vector, if needed.
        if ((maxUse == -1) || (parentESTidx == -1)) {
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
    // All done with no problems...
    return NO_ERROR;
}

int
MSTClusterMaker::worker() {
    // In the first step the worker first receives the index of the
    // EST was just added to the MST.  If the estNode added is -1,
    // then the MST building process is done.
    int estAdded = -1;
    // Wait for the Manager to send requests to this worker to perform
    // different tasks.
    MPI::Status msgInfo;
    do {
        // Wait for manager to send us a work request.  Since we are
        // waiting it should be tracked under idle time.
        MPI_PROBE(MANAGER_RANK, MPI_ANY_TAG, msgInfo);
        if (msgInfo.Get_tag() == COMPUTE_SIMILARITY_REQUEST) {
            // Read the actual message first.  Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            int estIdx = -1;
            MPI_RECV(&estIdx, 1, MPI::INT, MANAGER_RANK,
                     COMPUTE_SIMILARITY_REQUEST);
            // Perform the necessary operations.
            populateCache(estIdx);
        } else if (msgInfo.Get_tag() == COMPUTE_MAX_SIMILARITY_REQUEST) {
            // Read the actual message first. Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            int dummy = 0;
            MPI_RECV(&dummy, 1, MPI::INT, MANAGER_RANK,
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
            MPI_SEND(bestEntry, 4, MPI::INT, MANAGER_RANK,
                     MAX_SIMILARITY_RESPONSE);
        } else if (msgInfo.Get_tag() == ADD_EST) {
            // The manager has broad casted the next est to be added.
            MPI_RECV(&estAdded, 1, MPI::INT, MANAGER_RANK, ADD_EST);
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
                     MPI_INT, MANAGER_RANK, REPOPULATE_REQUEST);
        }
    } while (estAdded != -1);
    // Everything went on without a hitch.
    return NO_ERROR;
}

int
MSTClusterMaker::getOwnerProcess(const int estIdx) const {
    const int ESTsPerProcess = EST::getESTList().size() /
        MPI::COMM_WORLD.Get_size();
    const int ExtraESTs      = EST::getESTList().size() %
        MPI::COMM_WORLD.Get_size();
    const int ExtraESTsCutOff= (ExtraESTs * ESTsPerProcess) + ExtraESTs;
    
    // If the estIdx is less that the ExtraESTsCutOff then account for
    // the fact that each of these workers have one extra ESTs (as the
    // number of ESTs may not be evenly divisible by number of
    // processes).
    if (estIdx < ExtraESTsCutOff) {
        return estIdx / (ESTsPerProcess + 1);
    }
    // This est is in a process that does not have one extra...
    return (estIdx - ExtraESTs) / ESTsPerProcess;
}

void
MSTClusterMaker::populateCache(const int estIdx) {
    // First determine the list of ESTs that this process must deal
    // with using the helper method.
    int startESTidx, endESTidx;
    getOwnedESTidx(startESTidx, endESTidx);
    // Now compute similarity metric and store information in a SMList
    // data structure.
    SMList smList;
    analyzer->setReferenceEST(estIdx);
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if ((otherIdx == estIdx) || (cache->isESTinMST(otherIdx))) {
            // This EST entry can be ignored as this similarity metric
            // is not needed for MST construction.
            continue;
        }
        // Get similarity/distance metric.
        const float metric = analyzer->analyze(otherIdx);
        ASSERT ( metric >= 0 );
        // Obtain any alignemnt data that the analyzer may have.
        int alignmentData = 0;
        analyzer->getAlignmentData(alignmentData);
        // Add the information to metric list
        smList.push_back(CachedESTInfo(otherIdx, metric, alignmentData));
    }
    // Sort and prune the SMList to place ESTs with maximum similarity
    // at the top.
    cache->sortAndPrune(smList, cacheSize);
    // Now further process smList...
    const int ownerRank = getOwnerProcess(estIdx);
    if (ownerRank != MPI::COMM_WORLD.Get_rank()) {
        // This process is not the owner.  In this case, just send the
        // smList to the remote owner process.  There is some fudging
        // of data types going on here using the assumption that
        // sizeof(int) == sizeof(float).  Maybe can make this code
        // more MPI friendly.
        
        // First ensure there is at least one entry to work with.
        if (smList.empty()) {
            smList.push_back(CachedESTInfo(-1, -1.0f, -1));
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
    cache->mergeList(estIdx, smList, cacheSize);
    // Obtain similarity lists from all other processes (other than
    // ourselves).
    const int MyRank = MPI::COMM_WORLD.Get_rank();
    for(int pid = 0; (pid < MPI::COMM_WORLD.Get_size()); pid++) {
        if (pid == MyRank) {
            // Can't get a message from ourselves...
            continue;
        }
        // Choose actual rank depending on strict ordering scheme..
        const int rank = (strictOrder ? pid : MPI_ANY_SOURCE);
        MPI::Status msgInfo;
        // Wait to receive similarity list.  Since we are waiting it
        // should be logged as idle time for this process.
        MPI_PROBE(rank, SIMILARITY_LIST, msgInfo);
        // OK, we have a valid similarity list pending from some other
        // process. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI::CHAR);
        SMList remoteList(dataSize / sizeof(CachedESTInfo));
        // The following call is a kludge with MPI/STL data types
        // based on several language assumptions.  This part could be
        // cleaned up to be more portable later on.
        MPI_RECV(&remoteList[0], dataSize, MPI::CHAR,
                 msgInfo.Get_source(), SIMILARITY_LIST);
        // Merge the list we got from the remote process with our
        // local cache information if the list has a valid entry.
        if (hasValidSMEntry(remoteList)) {
            cache->mergeList(estIdx, remoteList, cacheSize);
        }
    }
    // Now finally let the manager know that the round of similarity
    // computation is all completed (assuming this process itself is
    // not the manager).
    if (MPI::COMM_WORLD.Get_rank() != MANAGER_RANK) {
        const int dummy = -1;
        MPI_SEND(&dummy, 1, MPI_INT, MANAGER_RANK, SIMILARITY_COMPUTATION_DONE);
    }
}

void
MSTClusterMaker::getOwnedESTidx(int& startIndex, int& endIndex) {
    const int ESTsPerProcess = EST::getESTList().size() /
        MPI::COMM_WORLD.Get_size();
    const int ExtraESTs      = EST::getESTList().size() %
        MPI::COMM_WORLD.Get_size();
    const int MyRank         = MPI::COMM_WORLD.Get_rank();
    
    // First figure out the starting and ending EST this processs is
    // responsible for further use.
    startIndex = MyRank * ESTsPerProcess;
    // Check for extra proceses as needed.
    if (MyRank <= ExtraESTs) {
        // The previous processes have one extra ESTs as number of
        // ESTs are not evenly divisible by number of processes.  So
        // account for this.
        startIndex = ((ESTsPerProcess + 1) * MyRank);
    } else {
        startIndex += ExtraESTs;
    }
    
    // Compute the last est index this process owns.
    endIndex = startIndex + ESTsPerProcess;
    if (MyRank < ExtraESTs) {
        // This process owns one extra EST as ESTs are not evenly
        // divisible by the number of processes.
        endIndex++;
    }
}

int
MSTClusterMaker::makeClusters() {
    int result = NO_ERROR;

    // All the parallel MPI processes load the data for processing
    // into memory.
    if ((result = analyzer->initialize()) != NO_ERROR) {
        // Error occured during initialization. Bail out.
        return result;
    }

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
        int startESTidx, endESTidx;
        getOwnedESTidx(startESTidx, endESTidx);
        cache = new MSTCache(EST::getESTList().size(), startESTidx,
                             endESTidx - startESTidx, analyzer,
                             !noCacheRepop);
        // Act as manager or worker depending on MPI rank.
        if (MPI::COMM_WORLD.Get_rank() == MANAGER_RANK) {
            // Get this MPI process to act as the manager.
            result = manager();
        } else {
            // Since rank in non-zero, this process must behave as a
            // worker
            result = worker();
        }
        // Dump cache usage statistics for this process.
        cache->displayStats(std::cout, MPI::COMM_WORLD.Get_rank());
        // Display MPI usage statistics.
        MPIStats::displayStats(std::cout);
        // Delete the cache as we no longer needed it.
        delete cache;
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
    // Do clustering if so desired.
    if (!dontCluster) {
        // Now get the helper to build the clusters.
        MSTCluster root;
        const double threshold = root.makeClusters(mst->getNodes(), percentile);
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
MSTClusterMaker::sendToWorkers(int data, const int tag) const {
    const int ProcessCount = MPI::COMM_WORLD.Get_size();
    for(int rank = 1; (rank < ProcessCount); rank++) {
        MPI_SEND(&data, 1, MPI_INT, rank, tag);
    }
}

bool
MSTClusterMaker::hasValidSMEntry(const SMList& list) const {
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
