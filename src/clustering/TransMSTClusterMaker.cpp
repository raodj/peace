#ifndef TRANS_MST_CLUSTER_MAKER_CPP
#define TRANS_MST_CLUSTER_MAKER_CPP

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

#include "TransMSTClusterMaker.h"
#include "HeuristicFactory.h"
#include "TVHeuristic.h"
#include "MPIStats.h"

TransMSTClusterMaker::TransMSTClusterMaker(ESTAnalyzer *analyzer) :
    MSTClusterMaker("tmst", analyzer) {
    // We don't have a reference est yet.
    currRefESTidx = -1;
    // Initialize counters to zero
    analyzeCount  = 0;
    successCount  = 0;
    // Initialize badMetric based on analyzer
    badMetric = analyzer->getInvalidMetric();
}

TransMSTClusterMaker::~TransMSTClusterMaker() {
    int startESTidx, endESTidx;
    getLocallyOwnedESTidx(startESTidx, endESTidx);
    for (int i = startESTidx; (i < endESTidx); i++) {
        if (metricCache[i] != NULL) {
            delete metricCache[i];
        }
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

bool
TransMSTClusterMaker::initialize() {
    if (!MSTClusterMaker::initialize()) {
        // Base class initialization failed. Major booboo.
        return false;
    }
    if (analyzer->getHeuristicChain() == NULL) {
        std::cerr << "TransMSTClusterMaker: No suitable heuristic found!"
                  << "\n";
        return false;
    }
    
    // Reserve space for TransCacheEntrys
    ASSERT ( estList != NULL );
    const size_t EstListSize = estList->size();
    metricCache.reserve(EstListSize);
    metricCache.insert(metricCache.begin(), EstListSize, NULL);
    // Everything went well.
    return true;
}

float
TransMSTClusterMaker::analyze(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    // Track number of calls to this method
    analyzeCount++;
    // Run analysis without using the heavy weight metric
    if (analyzer->compareMetrics(analyzer->analyze(otherEST, true, false),
                                 badMetric)) {
        ASSERT ( currRefESTidx >= 0 );
        // See if we have a main entry and a metric for otherEST
        TransCacheEntry* mainEntry = metricCache[otherEST->getID()];
        float transMetric;
        if (mainEntry != NULL &&
            mainEntry->getMetric(currRefESTidx, transMetric)) {
            // Yes! use metric from applying transitivity. Track
            // number of successes.
            successCount++;
            return transMetric;                
        }
        // When control drops here that means we have to use the heavy
        // weight analyzer for our operation.
        return analyzer->analyze(otherEST, false);
    }
    // When control drops here it indicates that heuristics failed.
    return badMetric;
}

void
TransMSTClusterMaker::pruneMetricEntries(const SMList& list,
                                         SMList& goodEntries) {
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
        CachedESTInfo dummy(-1, -1, badMetric, -1, -1);
        goodEntries.push_back(dummy);
    }
}

void
TransMSTClusterMaker::populateCache(const int estIdx,
                                   SMList* UNREFERENCED_PARAMETER(metricList)) {
    currRefESTidx = estIdx;
    // Now let the base class do its thing and populate the cache
    // entries for estIdx.
    SMList smList;
    MSTClusterMaker::populateCache(estIdx, &smList);
    // OK, now the cache entries for estIdx has been created at the
    // owner of this EST.  The owner process should send the list to
    // all the processes.
    const int ownerRank = getOwnerProcess(estIdx);
    if (ownerRank == MPI_GET_RANK()) {
        // This is the owner process. Send the cache entries to all
        // the processes.  First remove unwanted entries and get just
        // good entries.
        SMList goodEntries;
        pruneMetricEntries(smList, goodEntries);
        const int MsgSize      = goodEntries.size() * sizeof(CachedESTInfo);
        const int ProcessCount = MPI_GET_SIZE();

        const char *data       = reinterpret_cast<char*>(&goodEntries[0]);
        // An if check is necessary here when running on MPI
        // implementations that do not permit processes to
        // send/receive messages to themselves.
        for(int rank = 0; (rank < ProcessCount); rank++) {
            if (rank != ownerRank) {
                MPI_SEND(data, MsgSize, MPI_CHAR, rank, TRANSITIVITY_LIST);
            }
        }
        if ((goodEntries.size() == 1) && (goodEntries[0].estIdx == -1)) {
            // This is just an empty list. So don't process it.
            return;
        }
        processMetricList(goodEntries);
        // Finally, prune entries from our cache for the est that was
        // added in this method.
        delete metricCache[estIdx];
        metricCache[estIdx] = NULL;
    } else {
        // When control drops here, all the processes new receive the
        // SMList from the latest round of analysis. First probe to find
        // out size of the SMList.
        MPI_STATUS msgInfo;
        MPI_PROBE(ownerRank, TRANSITIVITY_LIST, msgInfo);
        // Now allocate the necessary space and obtain the SMList
        const int dataSize = MPI_GET_COUNT(msgInfo, MPI_TYPE_CHAR);
        SMList remoteList(dataSize / sizeof(CachedESTInfo));
        TRACK_IDLE_TIME(MPI_RECV(&remoteList[0], dataSize, MPI_CHAR,
                                 ownerRank, TRANSITIVITY_LIST));
        // Check to ensure that this not just a dummy list with one
        // invalid entry. Such dummy lists are required to work around MPI
        // implementations that do not permit zero-length messages.
        if ((remoteList.size() == 1) && (remoteList[0].estIdx == -1)) {
            // This is just an empty list. So don't process it.
            return;
        }
        processMetricList(remoteList);
    }
}

void
TransMSTClusterMaker::processMetricList(SMList& metricList) {
    // Now process the list of metrics we just received and build
    // caches for rapidly computing transitivity information. First
    // determine the set of ESTs this process owns
    int startESTidx, endESTidx;
    getLocallyOwnedESTidx(startESTidx, endESTidx);
    // Now process the list of metrics we just received and build
    // caches for rapidly computing transitivity information.
    for(size_t i = 0; (i < metricList.size()); i++) {
        // Add entries to our main metric cache
        if (metricList[i].estIdx >= startESTidx &&
            metricList[i].estIdx < endESTidx) {
            TransCacheEntry* entry = metricCache[metricList[i].estIdx];
            if (entry == NULL) {
                entry = new TransCacheEntry(metricList[i].estIdx);
            }
            entry->addEntries(metricList[i], metricList,
                              startESTidx, endESTidx);
            metricCache[metricList[i].estIdx] = entry;
        }
    }
}

#endif
