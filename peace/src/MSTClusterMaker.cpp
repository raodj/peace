#ifndef MST_CLUSTER_MAKER_CPP
#define MST_CLUSTER_MAKER_CPP

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

#include "MSTClusterMaker.h"
#include "ESTAnalyzer.h"
#include "MSTMultiListCache.h"
#include "MSTHeapCache.h"
#include "EST.h"
#include "MPIStats.h"
#include "HeuristicChain.h"
#include "Heuristic.h"
#include <fstream>
#include <sstream>

// A define to remove magic 0 (zero) in code
#define MANAGER_RANK 0

// Another define to remove magic 0 (zero) success checks.
#define NO_ERROR 0

// The default cache type to be used
char DefCacheType[10] = "heap";

// The static instance variables for command line arguments.
int    MSTClusterMaker::cacheSize     = 128;
bool   MSTClusterMaker::noCacheRepop  = true;
bool   MSTClusterMaker::strictOrder   = false;
char*  MSTClusterMaker::inputMSTFile  = NULL;
char*  MSTClusterMaker::outputMSTFile = NULL;
bool   MSTClusterMaker::dontCluster   = false;
bool   MSTClusterMaker::prettyPrint   = false;
bool   MSTClusterMaker::guiPrint      = false;
int    MSTClusterMaker::maxUse        = -1;
float  MSTClusterMaker::clsThreshold  = 1.0;
char*  MSTClusterMaker::cacheType     = DefCacheType;
char*  MSTClusterMaker::progFileName  = NULL;
    
// The common set of arguments for all FW EST analyzers
arg_parser::arg_record MSTClusterMaker::argsList[] = {
    {"--cache", "#similarity metrics to cache per EST",
     &MSTClusterMaker::cacheSize, arg_parser::INTEGER},
    {"--no-cache-repop", "Suppress EST cache repopulation",
     &MSTClusterMaker::noCacheRepop, arg_parser::BOOLEAN},    
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
    {"--gui-print", "Print the cluster tree for GUI processing.",
     &MSTClusterMaker::guiPrint, arg_parser::BOOLEAN},
    {"--maxUse", "Set a threshold to aggressively use metrics (default=0)",
     &MSTClusterMaker::maxUse, arg_parser::INTEGER},
    {"--clsThreshold", "Set a threshold for clustering (default=1.0)",
     &MSTClusterMaker::clsThreshold, arg_parser::FLOAT},
    {"--cacheType", "Set type of cache (heap or mlist) to use (default=heap)",
     &MSTClusterMaker::cacheType, arg_parser::STRING},
    {"--progress", "Log MST construction progress in a file (used by GUI)",
     &MSTClusterMaker::progFileName, arg_parser::STRING},        
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

MSTClusterMaker::MSTClusterMaker(ESTAnalyzer *analyzer,
                                 const int refESTidx,
                                 const std::string& outputFile)
    : ClusterMaker("mst", analyzer, refESTidx, outputFile), mst(NULL),
      hintKey_MST_RC("MST_RC") {
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
    const int WorkerCount = MPI_GET_SIZE();
    for(int workerID = 1; (workerID < WorkerCount); workerID++) {
        // Wait for a message to be received from a worker and obtain
        // some status information regarding the message.
        MPI_STATUS msgInfo;
        MPI_CODE({
                const int sourceRank=(strictOrder ? workerID : MPI_ANY_SOURCE);
                MPI_PROBE(sourceRank, REPOPULATE_REQUEST, msgInfo);
            });
        // OK, we have a valid repopulation request pending from some
        // worker. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI_TYPE_INT);
        int *requestData = new int[dataSize];
        MPI_RECV(requestData, dataSize, MPI_TYPE_INT,
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
            MPI_CODE({
                    int dummy;
                    TRACK_IDLE_TIME(MPI_RECV(&dummy, 1, MPI_TYPE_INT, ownerRank,
                                             SIMILARITY_COMPUTATION_DONE));
                });
        }
    }
    // Everything went on fine.
    return 0;
}

void
MSTClusterMaker::computeNextESTidx(int& parentESTidx, int& estToAdd,
                                   float& similarity, int& alignmentData,
                                   int& directionData) const {
    // First send request to compute local best EST values to each
    // worker process.
    sendToWorkers(-1, COMPUTE_MAX_SIMILARITY_REQUEST);
    // Now compute local best similarity
    cache->getBestEntry(parentESTidx, estToAdd, similarity, alignmentData,
                        directionData);
    // Receive similarity entry from 
    const int ProcessCount = MPI_GET_SIZE();
    for(int rank = 1; (rank < ProcessCount); rank++) {
        // Get the local simlarity information from another worker.
	int remoteData[5] = {0, 0, 0, 0, 0};
        MPI_CODE({
                // Choose worker rank depending on strict ordering scheme..
                const int workerRank = (strictOrder ? rank : MPI_ANY_SOURCE);
                TRACK_IDLE_TIME(MPI_RECV(remoteData, 5, MPI_TYPE_INT,
                                         workerRank, MAX_SIMILARITY_RESPONSE));
            });
        // Undo the fudge on similarity done at the sender end.
        const float remoteSim = *(reinterpret_cast<float*>(remoteData + 2));
        // Use a better or first valid entry
        if (analyzer->compareMetrics(remoteSim, similarity) ||
            (estToAdd == -1)) {
            // Found a higher similarity or a shorter distance in a
            // remote process or this is the first valid entry thusfar
            similarity    = remoteSim;
            parentESTidx  = remoteData[0];
            estToAdd      = remoteData[1];
            alignmentData = remoteData[3];
            directionData = remoteData[4];
        }
    }
}

void
MSTClusterMaker::addMoreChildESTs(const int parentESTidx, int& estToAdd,
                                  float &metric, int& alignmentData,
                                  int& directionData, int& pendingESTs) {
    // variable to track parent of next EST to be added.
    int newParent = -1;
    
    do {
        // Add the current node to the MST.
        ASSERT ( estToAdd != -1 );
        mst->addNode(parentESTidx, estToAdd, metric, alignmentData,
                     directionData);
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
        int childDirection = 0;
        computeNextESTidx(newParent, newChild, childMetric, childAlignment,
                          childDirection);
        // Check if we have a child with a useful alignment data.
        if ((newParent == parentESTidx) &&
            (analyzer->compareMetrics(childMetric, (float) maxUse))) {
            // Found a child of the same parent with a good metric. So
            // update parameters with the new child.
            estToAdd      = newChild;
            metric        = childMetric;
            alignmentData = childAlignment;
            directionData = childDirection;
        } else {
            // No more children to add!
            newParent = -1;
        }
    } while (newParent == parentESTidx);
}

void
MSTClusterMaker::updateProgress(const int estsAnalyzed,
                                const int totalESTcount) {
    if (progFileName == NULL) {
        // No need to report progress
        return;
    }
    // Open the progress file in a lazy manner as needed.
    if (!progressFile.is_open()) {
        progressFile.open(progFileName);
    }
    // Log the progress information.
    if (progressFile.good()) {
        progressFile << estsAnalyzed << "," << totalESTcount
                     << "\n" << std::flush;
        progressFile.seekp(0);
    }
}

int
MSTClusterMaker::manager() {
    // The number of pending nodes to be added to the MST.
    const int TotalESTcount = EST::getESTList().size() -
        EST::getProcessedESTCount();
    int pendingESTs         = TotalESTcount;
    // The minimum spanning tree that is built by this manager.
    int dummy;
    mst = new MST(pendingESTs, analyzer->getAlignmentData(dummy));
    // Kick off all activities by adding the reference EST as the root
    // of the MST with a similarity metric of 0.
    int parentESTidx    = -1;
    int estToAdd        = refESTidx;
    // Ensure the refESTidx has not been proceed out. If it has been
    // try to make best effort to find an alternative.
    if (EST::getEST(estToAdd)->hasBeenProcessed()) {
        std::cerr << "Warning: The reference EST has been filtered out.\n"
                  << "Trying to select the first non-filtered EST instead.\n";
        estToAdd = 0;
        while ((estToAdd < EST::getESTCount()) &&
               (EST::getEST(estToAdd)->hasBeenProcessed())) {
            estToAdd++;
        }
        if (estToAdd >= EST::getESTCount()) {
            std::cerr << "Error: All the ESTs have been filtered out.\n"
                      << "Nothing left to be clustered.\n";
            return 0;
        } else {
            std::cerr << "Using EST at index " << estToAdd
                      << " as root.\n";
        }
    }
    
    float metric        = analyzer->getValidMetric();
    int   alignmentInfo = 0;
    int   directionInfo = 0;
    do {
        // Update progress information as needed
        updateProgress(TotalESTcount - pendingESTs, TotalESTcount);
        // Add the EST to the MST vector, if needed.
        if ((maxUse == -1) || (parentESTidx == -1)) {
            ASSERT ( estToAdd != -1 );
            mst->addNode(parentESTidx, estToAdd, metric, alignmentInfo,
                         directionInfo);
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
        computeNextESTidx(parentESTidx, estToAdd, metric, alignmentInfo,
                          directionInfo);
        ASSERT( parentESTidx != -1 );
        ASSERT( estToAdd     != -1 );
        if (maxUse != -1) {
            // Try to add as many ESTs as possible rooted at the given
            // parentESTidx using a helper method.
            addMoreChildESTs(parentESTidx, estToAdd, metric,
                             alignmentInfo, directionInfo, pendingESTs);
        }
    } while (pendingESTs > 0);

    // Update progress information as needed
    updateProgress(TotalESTcount - pendingESTs, TotalESTcount);
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
    MPI_STATUS msgInfo;
    do {
        // Wait for manager to send us a work request.  Since we are
        // waiting it should be tracked under idle time.
        MPI_PROBE(MANAGER_RANK, MPI_ANY_TAG, msgInfo);
        if (msgInfo.Get_tag() == COMPUTE_SIMILARITY_REQUEST) {
            // Read the actual message first.  Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            int estIdx = -1;
            MPI_RECV(&estIdx, 1, MPI_TYPE_INT, MANAGER_RANK,
                     COMPUTE_SIMILARITY_REQUEST);
            // Perform the necessary operations.
            populateCache(estIdx);
        } else if (msgInfo.Get_tag() == COMPUTE_MAX_SIMILARITY_REQUEST) {
            // Read the actual message first. Dont' account it under
            // idle time as we are doing this Recv because the message
            // has already arrived.
            MPI_CODE({
                    int dummy = 0;
                    MPI_RECV(&dummy, 1, MPI_TYPE_INT, MANAGER_RANK,
                             COMPUTE_MAX_SIMILARITY_REQUEST);
                });
            int   bestEntry[5];
            float similarity = 0;
            // Get the best possible local similarity match.
            cache->getBestEntry(bestEntry[0], bestEntry[1],
                                similarity, bestEntry[3], bestEntry[4]);
            // Fudge the similarity into the bestEntry array to
            // transmitt the necessary information to the manager.
            // Maybe there is a cleaner way to do it too...
            int *temp    = reinterpret_cast<int*>(&similarity);
            bestEntry[2] = *temp;
            MPI_SEND(bestEntry, 5, MPI_TYPE_INT, MANAGER_RANK,
                     MAX_SIMILARITY_RESPONSE);
        } else if (msgInfo.Get_tag() == ADD_EST) {
            // The manager has broad casted the next est to be added.
            MPI_RECV(&estAdded, 1, MPI_TYPE_INT, MANAGER_RANK, ADD_EST);
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
                     MPI_TYPE_INT, MANAGER_RANK, REPOPULATE_REQUEST);
        }
    } while (estAdded != -1);
    // Everything went on without a hitch.
    return NO_ERROR;
}

int
MSTClusterMaker::getOwnerProcess(const int estIdx) const {
    const int ESTsPerProcess = EST::getESTList().size() / MPI_GET_SIZE();
    const int ExtraESTs      = EST::getESTList().size() % MPI_GET_SIZE();
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

float
MSTClusterMaker::analyze(const int otherEST) {
    return analyzer->analyze(otherEST);
}

int
MSTClusterMaker::initialize() {
    // All the parallel MPI processes load the data for processing
    // into memory.
    int result;
    if ((result = analyzer->initialize()) != NO_ERROR) {
        // Error occured during initialization. Bail out.
        return result;
    }
    // Return success.
    return NO_ERROR;
}

void
MSTClusterMaker::populateCache(const int estIdx, SMList* metricList) {
    // First determine the list of ESTs that this process must deal
    // with using the helper method.
    int startESTidx, endESTidx;
    getOwnedESTidx(startESTidx, endESTidx);
    // Setup the reference estIdx in the analyzer which given the
    // analyzer a chance to optimize initialization.
    analyzer->setReferenceEST(estIdx);
    // Now compute similarity metric and store information in a SMList
    // data structure.
    SMList smList;
    const float InvalidMetric = analyzer->getInvalidMetric();
    // Flag to ensure only one invalid metric gets added
    bool needInvalidMetric = true;
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if (cache->isESTinMST(otherIdx) || (otherIdx == estIdx) ||
            (EST::getEST(otherIdx)->hasBeenProcessed())) {
            // This EST entry can be ignored as this metric is not
            // needed (or must not be used) for MST construction.
            continue;
        }
        // Get similarity/distance metric.
        const float metric = analyze(otherIdx);
        ASSERT ( metric >= 0 );
        // Obtain any alignment data that the analyzer may have.
        int alignmentData = 0;
        int directionData = 1;
        analyzer->getAlignmentData(alignmentData);
        HeuristicChain* chain = HeuristicChain::getHeuristicChain();
        if (chain != NULL) {
            chain->getHint(HeuristicChain::MST_RC, directionData);
        }
        // Add only the first invalid entry. One is enough to do build
        // a valid MST. There is no need for multiple vestigial
        // entries.
        const bool isMetricOK= analyzer->compareMetrics(metric, InvalidMetric);
        if ((isMetricOK) || (needInvalidMetric)) {
            // Add the information to metric list
            smList.push_back(CachedESTInfo(estIdx, otherIdx,
                                           metric, alignmentData,
                                           directionData));
            // If this is an invalid entry, then one is enough. Don't
            // add more.
            needInvalidMetric= (needInvalidMetric && isMetricOK);
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
            smList.push_back(CachedESTInfo(-1, -1, -1.0f, -1, -1));
        }
        MPI_SEND(&smList[0], smList.size() * sizeof(CachedESTInfo),
                 MPI_TYPE_CHAR, ownerRank, SIMILARITY_LIST);
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
    const int MyRank = MPI_GET_RANK();
    for(int pid = 0; (pid < MPI_GET_SIZE()); pid++) {
        if (pid == MyRank) {
            // Can't get a message from ourselves...
            continue;
        }
        // The status structure for probing incoming messages.
        MPI_STATUS msgInfo;
        MPI_CODE({
                // Choose actual rank depending on strict ordering scheme..
                const int rank = (strictOrder ? pid : MPI_ANY_SOURCE);
                // Wait to receive similarity list.  Since we are waiting it
                // should be logged as idle time for this process.
                MPI_PROBE(rank, SIMILARITY_LIST, msgInfo);
            });
        // OK, we have a valid similarity list pending from some other
        // process. So read and process it.
        const int dataSize = msgInfo.Get_count(MPI_TYPE_CHAR);
        SMList remoteList(dataSize / sizeof(CachedESTInfo));
        // The following call is a kludge with MPI/STL data types
        // based on several language assumptions.  This part could be
        // cleaned up to be more portable later on.
        MPI_RECV(&remoteList[0], dataSize, MPI_TYPE_CHAR,
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
    if (MPI_GET_RANK() != MANAGER_RANK) {
        MPI_CODE({
                const int dummy = -1;
                MPI_SEND(&dummy, 1, MPI_TYPE_INT, MANAGER_RANK,
                         SIMILARITY_COMPUTATION_DONE);
            });
    }
}

void
MSTClusterMaker::getOwnedESTidx(int& startIndex, int& endIndex) {
    const int ESTsPerProcess = EST::getESTList().size() / MPI_GET_SIZE();
    const int ExtraESTs      = EST::getESTList().size() % MPI_GET_SIZE();
    const int MyRank         = MPI_GET_RANK();
    
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

void
MSTClusterMaker::displayStats(std::ostream& os) {
    // Dump cache usage statistics for this process.
    cache->displayStats(os, MPI_GET_RANK());
    // Display MPI usage statistics.
    MPIStats::displayStats(os);
}

int
MSTClusterMaker::populateMST() {
    if (inputMSTFile != NULL) {
        // In this case the MST is to be read from an input file.
        // So let's do that.
        if ((mst = MST::deSerialize(inputMSTFile)) == NULL) {
            // Some error occured loading MST. Can't proceed further
            return 2;
        }
    }
    
    // No input MST file supplied.  That means the MST must be
    // built. The local cache that contains information to build the
    // MST.  The cache for that contains information regarding ESTs
    // under the ownership of this process.  The cache facilitates
    // rapid construction of ESTs but minimizing the number of times
    // similarity metrics need to computed.  The cache needs to know
    // the starting indedx of the EST whose data must be cached by
    // this process.  So this information is computed first.
    int result = NO_ERROR;
    int startESTidx, endESTidx;
    getOwnedESTidx(startESTidx, endESTidx);
    // Create a suitable cache.
    if (!strcmp(cacheType, "mlist")) {
        cache = new MSTMultiListCache(EST::getESTList().size(), startESTidx,
                                      endESTidx - startESTidx, analyzer,
                                      !noCacheRepop, cacheSize);
    } else {
        cache = new MSTHeapCache(EST::getESTList().size(), startESTidx,
                                 endESTidx - startESTidx, analyzer,
                                 !noCacheRepop, cacheSize);
    }
    ASSERT ( cache != NULL );
    // Act as manager or worker depending on MPI rank.
    if (MPI_GET_RANK() == MANAGER_RANK) {
        // Get this MPI process to act as the manager.
        result = manager();
    } else {
        // Since rank in non-zero, this process must behave as a
        // worker
        result = worker();
    }
    // Now either we have a MST built or there is an error.
    return result;
}

int
MSTClusterMaker::addDummyCluster(const std::string name) {
    MSTCluster *dummy = new MSTCluster(&root, name);
    root.add(dummy);
    // Return id of newly added cluster
    return dummy->getClusterID();
}

void
MSTClusterMaker::addEST(const int clusterID, const int estIdx) {
    // Create a dummy MSTNode.
    MSTNode node(-1, estIdx, analyzer->getInvalidMetric());
    root.add(clusterID, node);
    // Flag the EST as having been processed
    EST::getEST(estIdx)->setProcessed(true);
}

int
MSTClusterMaker::buildAndShowClusters() {
    // Let the root cluster build sub-clusters
    root.makeClusters(mst->getNodes(), analyzer, clsThreshold);
    
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
    // Choose output stream depending on condition of outFile.
    std::ostream& dest = (outFile.is_open() && outFile.good()) ?
        outFile : std::cout;
    if (prettyPrint) {
        root.printClusterTree(dest);
    } else if (guiPrint) {
        root.guiPrintClusterTree(dest, analyzer->getInputFileName());
    }  else {
        // No pretty printing. Just dump the info out.
        dest << root << std::endl;
    }
    // Can't have errors (yet in this method).
    return NO_ERROR;
}


int
MSTClusterMaker::makeClusters() {
    int result = NO_ERROR;

    // Next compute/load MST using helper method.
    result = populateMST();
    
    // Display statistics by writing all the data to a string stream
    // and then flushing the stream.  This tries to working around (it
    // is not perfect solution) interspersed data from multiple MPI
    // processes
    std::ostringstream buffer;
    displayStats(buffer);
    std::cout << buffer.str() << std::endl;
    // Delete the cache as we no longer needed it.
    delete cache;

    if ((result != NO_ERROR) || (mst == NULL)) {
        // Some error occured.  Can't proceed further.
        return result;
    }
    
    ASSERT ( mst != NULL );
    // Dump MST out to the specified output file (if any).
    if ((outputMSTFile != NULL) && (mst != NULL)) {
        mst->serialize(outputMSTFile, (inputMSTFile != NULL) ? inputMSTFile :
                       analyzer->getInputFileName(), clsThreshold);
    }
    
    // Do clustering and display results if so desired.
    if (!dontCluster) {
        result = buildAndShowClusters();
    }
    // Return result from behaving as a manager or a worker.
    return result;
}

void
MSTClusterMaker::sendToWorkers(int data, const int tag) const {
    const int ProcessCount = MPI_GET_SIZE();
    for(int rank = 1; (rank < ProcessCount); rank++) {
        MPI_SEND(&data, 1, MPI_TYPE_INT, rank, tag);
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
