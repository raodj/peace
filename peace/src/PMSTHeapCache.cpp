#ifndef PMST_HEAP_CACHE_CPP
#define PMST_HEAP_CACHE_CPP

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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "PMSTHeapCache.h"

PMSTHeapCache::PMSTHeapCache(const int totalESTCount,
                             const int startESTidx,
                             const int numOwnESTs,
                             const ESTAnalyzer *estAnalyzer):
    MSTHeapCache(totalESTCount, startESTidx, numOwnESTs, estAnalyzer,
                 false, totalESTCount) {
    // Nothing else to be done for now in the constructor.
}

void
PMSTHeapCache::popBestEntry(int& srcESTidx, int& destESTidx,
                            float& metric, int& alignmentData) {
    getBestEntry(srcESTidx, destESTidx, metric, alignmentData);
    // Make sure we found an entry
    if (srcESTidx != -1) {
        // Remove that entry from the cache
        popCache();
    }
}

#endif
