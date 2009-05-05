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
#include <iterator>

MSTCache::MSTCache(const int totalESTCount, const int startESTidx,
                   const int numberOfOwnESTs,
                   const ESTAnalyzer *estAnalyzer,
                   const bool repopulate, const int maxCacheSize) :
    analyzer(estAnalyzer), nodesAdded(totalESTCount, false),
    startOwnESTidx(startESTidx), numOwnESTs(numberOfOwnESTs),
    maxCachePerEST(maxCacheSize), repopulateCaches(repopulate) {
    // Nothing else to be done for now in the constructor.
}

// static utility method
void
MSTCache::copy_n(const SMList& input, const size_t count, SMList &output) {
    // Get minimum of input.size() or count to ensure we never exceed
    // the number of elements in input list.
    const size_t stopIndex = (count < input.size()) ? count : input.size();
    // Copy the necessary elements from input to output...
    for(size_t index = 0; (index < stopIndex); index++) {
        output.push_back(input[index]);
    }
}

#endif
