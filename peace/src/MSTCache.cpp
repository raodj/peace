#ifndef MST_CACHE_CPP
#define MST_CACHE_CPP

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

#include "MSTCache.h"
#include "Utilities.h"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <cmath>
#include <mpi.h>

MSTCache::MSTCache(const int totalESTCount, const int startESTidx,
                   const int numOwnESTs,
                   const ESTAnalyzer *estAnalyzer) :
    analyzer(estAnalyzer), cache(numOwnESTs),
    nodesAdded(totalESTCount, false),
    startOwnESTidx(startESTidx), cacheRepopulation(0),
    prunedEntries(0), analysisCount(0) {
    // Nothing else to be done for now in the constructor.
}

void
MSTCache::pruneCaches(const int estIdxAdded,
                      std::vector<int>& repopulateList,
                      const bool prefixListSize) {
    // Clear out any entries in the list of ESTs whose cache
    // information needs to be recomputed.
    repopulateList.clear();
    if (prefixListSize) {
        // Reserve space for number of entries in the list.
        int temp = 0;
        repopulateList.push_back(temp);
    }
    
    // Prune each SMList in the cache using a helper method.
    for(std::vector<MSTCacheEntry>::iterator curr = cache.begin();
        (curr != cache.end()); curr++) {
        MSTCacheEntry& entry = *curr;
        if (entry.second.empty()) {
            // Nothing to prune if the list is already empty.
            continue;
        }
        if (pruneCache(entry.second, estIdxAdded)) {
            // If pruning the SMList made the list empty then it needs
            // to be repopulated.
            repopulateList.push_back(entry.first);
            // Track number of times cache had to be repouplated
            cacheRepopulation++;
        }
    }
    // Now, set the flag for the newly added est node.
    nodesAdded[estIdxAdded] = true;

    // Finally update the number of entries in the list as needed.
    if (prefixListSize) {
        // Update number of entries actually added to the list
        // ignoring the the size entry itself.
        repopulateList[0] = repopulateList.size() - 1;
    }
}

bool
MSTCache::pruneCache(SMList& list, const int estIdx) {
    SMList::iterator curr = list.begin();
    while (curr != list.end()) {
        if ((*curr).first == estIdx) {
            // This entry needs to be erased.  Track the number of
            // pruned entries
            prunedEntries++;
            // Erase entry and reset the iterator
            curr = list.erase(curr);
        } else {
            // On to the next SMEntry
            curr++;
        }
    }
    
    return list.empty();
}

void
MSTCache::mergeList(const int estIdx, const SMList& list,
                    const int maxCacheSize) {
    // If the list to be merged is empty there is absolutely nothing
    // left for this method to do.
    if (list.empty()) {
        // Nothing to be done.
        return;
    } 
    // Extract the necessary information to streamline the code
    const int localIndex      = estIdx - startOwnESTidx;
    MSTCacheEntry& cacheEntry = cache[localIndex];
    SMList& cacheList         = cacheEntry.second;
    
    // First compute a constant worst metric value depending on
    // whether similarity or distance metrics are being used. In the
    // case of similarity, worstMetric will be -1.
    float worstMetric = analyzer->getInvalidMetric();
    if (cacheList.size() > 0) {
        // Cache has entries.  Obtain the least useful (that is, most
        // dissimilar or most distant) entry in the cache.
        worstMetric = cacheList.back().second;
    } else {
        // Update the estIdx for this cache entry for the first time.
        cacheEntry.first = estIdx;
    }

    // First check if the merge operation needs to be done by checking
    // similarity/distance metrics of the two lists.
    SMList::const_iterator mergeStart = list.begin();
    if (analyzer->compareMetrics(worstMetric, (*mergeStart).second) &&
        ((int) cacheList.size() >= maxCacheSize)) {
        // All the entries in the cache have higher similarity metric
        // than the highest similarity in this list and the cache is
        // full.  So, nothing further to be done!
        return;
    }

    // Merge the two lists into a new list
    SMList mergedList(list.size() + cacheList.size());
    LessSMEntry mergeComparator(analyzer);
    std::merge(list.begin(), list.end(),
               cacheList.begin(), cacheList.end(),
               mergedList.begin(), mergeComparator);
    // Copy the appropriate number of entries from the merged list
    // back to the original list of entries.
    cacheList.clear();
    if ((int) mergedList.size() > maxCacheSize) {
        if ((int) cacheList.size() != maxCacheSize) {
            // Ensure cacheList is of the right size first.
            cacheList = SMList(maxCacheSize);
        }
        // Copy only the maxCacheSize entries...
        std::copy(mergedList.begin(), mergedList.begin() + maxCacheSize,
                  cacheList.begin());
    } else  {
        // Copy the whole merged set of entries
        cacheList = mergedList;
    }
    
    ASSERT ( cacheList.size() >= list.size() );
}

void
MSTCache::getBestEntry(int& srcESTidx, int& destESTidx,
                       float& metric) const {
    // Initialize the results to an invalid value first.
    srcESTidx= destESTidx = -1;
    metric   = analyzer->getInvalidMetric();
    
    // Locate the best entry from various lists in this cache.
    for(std::vector<MSTCacheEntry>::const_iterator curr = cache.begin();
        (curr != cache.end()); curr++) {
        const MSTCacheEntry& cacheEntry = *curr;
        if (cacheEntry.second.empty()) {
            // This list is empty.  So there is nothing to compare.
            continue;
        }
        // Get the entry with best metric (namely the first entry as
        // these lists are already sorted) in this list.
        const SMEntry& smEntry = cacheEntry.second.front();
        // Is this entry better than the ones before?
        if (analyzer->compareMetrics(smEntry.second, metric)) {
            // Yes, it is!  Good. update the necessary information.
            srcESTidx  = cacheEntry.first;
            destESTidx = smEntry.first;
            metric     = smEntry.second;
        }
    }
}

void
MSTCache::sortAndPrune(SMList& list, const int cacheSize) {
    // Sort the list such that the most similar entries are to the top
    // of the list so ease pruning less significant entries.
    std::sort(list.begin(), list.end(), LessSMEntry(analyzer));
    // Use number of entries in the list to determine number of
    // analyses performed
    analysisCount += list.size();
    // Prune the list to ensure that the number of entries don't
    // exceed the desired cache sizes.
    while ((int) list.size() > cacheSize) {
        // Prune the trailing list in smList as they now have the
        // least similar EST entries there.
        list.pop_back();
        // Track number of prunes
        prunedEntries++;
    }
}

void
MSTCache::displayStats(std::ostream &os, const int MyRank) const {
    // First count total number of MST entries in the cache.
    int total = 0;
    for(std::vector<MSTCacheEntry>::const_iterator curr = cache.begin();
        (curr != cache.end()); curr++) {
        total += (*curr).second.size();
    }
    // Note that when control drops here, prunedEntries is also
    // accounting for entries that were removed and added to the MST.
    // This is good use of the cache and must not be accounted for in
    // the prunedEntries.
    int adjustedPrunedEntries = prunedEntries - cache.size();
    if (adjustedPrunedEntries < 0) { 
        adjustedPrunedEntries = 0;
    }
    
    os << "Statistics on process with Rank " << MyRank << "\n"
       << "-------------------------------------------\n"
       << "    #ESTs managed by cache   : " << cache.size() << "\n"
       << "    Total remaining entries  : " << total << "\n"
       << "    #Pruned entries          : " << adjustedPrunedEntries
       << std::endl;
    // Compute mean 
    const double mean = (double) total / cache.size();
    // Now compute variation.
    double variance = 0;
    for(std::vector<MSTCacheEntry>::const_iterator curr = cache.begin();
        (curr != cache.end()); curr++) {
        const double deviation = mean - (*curr).second.size();
        variance += (deviation * deviation);
    }
    // Now compute standard deviation
    const double deviation = sqrt(variance);
    // Display the mean and standard deviation.
    os << "    Mean (Standard deviation): " << mean
       << " (" << deviation << ")\n";
    // Display number of times this cache had to be repopulated.
    os << "    Cache Repopulation count : " << cacheRepopulation << "\n"
       << "    #EST analyses performed  : " << analysisCount
       << std::endl;
}

#endif