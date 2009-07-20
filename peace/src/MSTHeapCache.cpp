#ifndef MST_HEAP_CACHE_CPP
#define MST_HEAP_CACHE_CPP

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

#include "MSTHeapCache.h"

MSTHeapCache::MSTHeapCache(const int totalESTCount,
                           const int startESTidx,
                           const int numOwnESTs,
                           const ESTAnalyzer *estAnalyzer,
                           const bool repopulate,
                           const int  maxCachePerEST) :
    MSTCache(totalESTCount, startESTidx, numOwnESTs, estAnalyzer,
             repopulate, maxCachePerEST), 
    prunedEntries(0), cache(GreaterCachedESTInfo(estAnalyzer)) {
    // Nothing else to be done for now in the constructor.
}

void
MSTHeapCache::pruneCaches(const int estIdxAdded,
                          std::vector<int>& repopulateList,
                          const bool prefixListSize) {
    // Clear out any entries in the list of ESTs whose cache
    // information needs to be recomputed.
    repopulateList.clear();
    if (prefixListSize) {
        // Reserve space for number of entries in the list.
        repopulateList.push_back(0);
    }
    // Now, set the flag for the newly added est node.
    nodesAdded[estIdxAdded] = true;

    // Clear out entries that correspond to the estIdxAdded from the
    // heap.
    while (!cache.empty() && (isESTinMST(cache.top().estIdx))) {
        popCache();
    }
}

void
MSTHeapCache::mergeList(const int UNREFERENCED_PARAMETER(estIdx), 
                        const SMList& list) {
    // Add all entries in the list to our cache.
    for(size_t i = 0; (i < list.size()); i++) {
        cache.push(list[i]);
    }
}

void
MSTHeapCache::getBestEntry(int& srcESTidx, int& destESTidx,
                           float& metric, int & alignmentData) const {
    // Initialize the results to an invalid value first.
    srcESTidx= destESTidx = alignmentData = -1;
    metric   = analyzer->getInvalidMetric();

    if (!cache.empty()) {
        // Simply use the top-element as the best choice
        const CachedESTInfo& entry = cache.top();
        srcESTidx     = entry.refESTidx;
        destESTidx    = entry.estIdx;
        metric        = entry.metric;
        alignmentData = entry.alignmentData;
    }
}

void
MSTHeapCache::displayStats(std::ostream &os, const int MyRank) const {
    const int totalEntries = cache.size() + prunedEntries;
    
    os << "Statistics on process with Rank " << MyRank << "\n"
       << "-------------------------------------------\n"
       << "    #ESTs entries managed    : " << totalEntries << "\n"
       << "    Total remaining entries  : " << cache.size() << "\n"
       << "    #Pruned entries          : " << prunedEntries
       << std::endl;
}

void
MSTHeapCache::popCache() {
    // Remove entry from top of the cache.
    cache.pop();
    // Track number of entries pruned
    prunedEntries++;
}

#endif
