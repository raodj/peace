#ifndef TRANS_MST_CLUSTER_MAKER_CPP
#define TRANS_MST_CLUSTER_MAKER_CPP

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

#include "TransMSTClusterMaker.h"
#include "HeuristicFactory.h"
#include "TVHeuristic.h"
#include "MPIStats.h"

TransMSTClusterMaker::TransMSTClusterMaker(ESTAnalyzer *analyzer,
                                           const int refESTidx,
                                           const std::string& outputFile) :
    MSTClusterMaker(analyzer, refESTidx, outputFile) {
    // Initialize tvHeuristic pointer that will be set later on.
    tvHeuristic   = NULL;
    // We don't have a reference est yet.
    currRefESTidx = 0;
    // Initialize counters to zero
    analyzeCount  = 0;
    successCount  = 0;
}

TransMSTClusterMaker::~TransMSTClusterMaker() {
    if (tvHeuristic != NULL) {
        // We have a valid tvHeuristic object. This object could have
        // been from the heuristicChain in which case, it must not be
        // deleted here by this class.
        if ((analyzer->getHeuristicChain() == NULL) ||
            (analyzer->getHeuristicChain()->getHeuristic("tv")!= tvHeuristic)) {
            // The tvHeuristic object was custom created and is not
            // from the heuristic chain.
            delete tvHeuristic;
        }
        // Reset pointer to detect any lurking memory errors
        tvHeuristic = NULL;
    }
}

void
TransMSTClusterMaker::displayStats(std::ostream& os) {
    // First let the base class display statistics.
    MSTClusterMaker::displayStats(os);
    // Now print additional stats on conditional-transitivity usage
    os << "    Total calls to analyze   : " << analyzeCount << "\n"
       << "    Successful transitivity  : " << successCount << "\n";
}

void
TransMSTClusterMaker::setTVHeuristic() {
    if (tvHeuristic != NULL) {
        // The pointer has already been set. Nothing further to do.
        return;
    }
    // First see if the analyzer's heuristic chain already has a tv
    // heuristic class. If so use it.
    if (analyzer->getHeuristicChain() != NULL) {
        tvHeuristic = analyzer->getHeuristicChain()->getHeuristic("tv");
        // Insanity check to ensure returned pointer is a tv heuristic
        // or a subclass of it.
        ASSERT ( dynamic_cast<TVHeuristic*>(tvHeuristic) != NULL );
    }
    // If the analyzer chain did not have a tv heuristic object, then
    // create one for our own use.
    if (tvHeuristic == NULL) {
        tvHeuristic = HeuristicFactory::create("tv", 0, "");
        ASSERT ( tvHeuristic != NULL );
    }
}

float
TransMSTClusterMaker::analyze(const int otherEST) {
    // Track number of calls to this method
    analyzeCount++;
    
    // First, let's see if our metricCache has an entry for other EST.
    ASSERT ( currRefESTidx >= 0 );
    // See if we have a main entry for otherEST
    MetricCacheMap::iterator mainEntry = metricCache.find(otherEST);
    if (mainEntry != metricCache.end()) {
        float transMetric;
        // See if the TransCacheEntry has a metric for otherEST
        if (mainEntry->second.getMetric(currRefESTidx, transMetric)) {
            // Yes! It does. OK. Do these two pass t/v heuristic?
            ASSERT ( tvHeuristic != NULL );
            if (tvHeuristic->shouldAnalyze(otherEST)) {
                // Yes! use metric from applying transitivity. Track
                // number of successes.
                successCount++;
                return transMetric;
            }
        }
    }
    // When control drops here that means we have to use the heavy
    // weight analyzer for our operation.
    return analyzer->analyze(otherEST);
}

void
TransMSTClusterMaker::pruneMetricEntries(const SMList& list,
                                         SMList& goodEntries,
                                         const float badMetric) {
    // Copy all the necessary entries to a temporary list as this make
    // the process faster.
    goodEntries.reserve(list.size());    
    for(SMList::const_iterator entry = list.begin();
        (entry != list.end()); entry++) {
        if ((analyzer->compareMetrics(entry->metric, badMetric))) {
            // This is a good entry. Save it.
            goodEntries.push_back(*entry);
        }
    }
    // If the list has zero entries, then add at least one invalid
    // entry to ensure that the list is not empty.
    if (goodEntries.empty()) {
        CachedESTInfo dummy(-1, -1, badMetric, -1);
        goodEntries.push_back(dummy);
    }
}

void
TransMSTClusterMaker::populateCache(const int estIdx,
                                    SMList* UNREFERENCED_PARAMETER(list)) {
    // First setup our current reference EST with the heuristic for
    // its use and setup.
    setTVHeuristic();
    ASSERT ( tvHeuristic != NULL );
    tvHeuristic->setReferenceEST((currRefESTidx = estIdx));
    // Now let the base class do its thing and populate the cache
    // entries for estIdx.
    SMList smList;
    MSTClusterMaker::populateCache(estIdx, &smList);
    // OK, now the cache entries for estIdx has been created at the
    // owner of this EST.  The owner process should send the list to
    // all the processes.
    const int ownerRank = getOwnerProcess(estIdx);
    if (ownerRank == MPI::COMM_WORLD.Get_rank()) {
        // This is the owner process. Send the cache entries to all
        // the processes.  First remove unwanted entries and get just
        // good entries.
        SMList goodEntries;
        pruneMetricEntries(smList, goodEntries, analyzer->getInvalidMetric());
        const int MsgSize      = goodEntries.size() * sizeof(CachedESTInfo);
        const int ProcessCount = MPI::COMM_WORLD.Get_size();
        const char *data       = reinterpret_cast<char*>(&goodEntries[0]);
        for(int rank = 0; (rank < ProcessCount); rank++) {
            MPI_SEND(data, MsgSize, MPI_CHAR, rank, TRANSITIVITY_LIST);
        }
    }
    
    // When control drops here, all the processes new receive the
    // SMList from the latest round of analysis. First probe to find
    // out size of the SMList.
    MPI::Status msgInfo;
    MPI_PROBE(ownerRank, TRANSITIVITY_LIST, msgInfo);
    // Now allocate the necessary space and obtain the SMList
    const int dataSize = msgInfo.Get_count(MPI::CHAR);
    if (dataSize <= 0) {
        // There are no remote cache entries to process.
        return;
    }
    SMList remoteList(dataSize / sizeof(CachedESTInfo));
    MPI_RECV(&remoteList[0], dataSize, MPI::CHAR, ownerRank,TRANSITIVITY_LIST);
    // Check to ensure that this not just a dummy list with one
    // invalid entry. Such dummy lists are required to work around MPI
    // implementations that do not permit zero-length messages.
    if ((remoteList.size() == 1) && (remoteList[0].estIdx == -1)) {
        // This is just an empty list. So don't process it.
        return;
    }
    // Now process the list of metrics we just received and build
    // caches for rapidly computing transitivity information. First
    // determine the set of ESTs this process owns
    int startESTidx, endESTidx;
    getOwnedESTidx(startESTidx, endESTidx);
    for(size_t i = 0; (i < remoteList.size()); i++) {
        // Add entries to our main metric cache
        TransCacheEntry &entry = metricCache[remoteList[i].estIdx];
        entry.estIdx           = remoteList[i].estIdx;
        entry.addEntries(remoteList[i], remoteList, startESTidx, endESTidx);
    }
    // Finally, prune entries from our cache for the est that was
    // added in this method.
    metricCache.erase(estIdx);
}

#endif
