#ifndef TRANS_CACHE_ENTRY_CPP
#define TRANS_CACHE_ENTRY_CPP

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

#include "TransCacheEntry.h"
#include "Utilities.h"

TransCacheEntry::TransCacheEntry(const int refEstIdx) : estIdx(refEstIdx) {
    peerInfo = new TransCacheMap;
    // Currently we have nothing else to be done here.
}

TransCacheEntry::~TransCacheEntry() {
    delete peerInfo;
}

void
TransCacheEntry::addEntries(const CachedESTInfo& reference,
                            const SMList& metrics,
                            const int UNREFERENCED_PARAMETER(startIndex), 
                            const int UNREFERENCED_PARAMETER(endIndex)) {
    // Add the related metric information to the peerInfo map
    for(size_t index = 0; (index < metrics.size()); index++) {
        const CachedESTInfo& entry = metrics[index];
        // Don't add entries referring to ourselves.
        if (entry.estIdx != estIdx) {
            // Metric is minimum of peer metric and parent metric,
            // as according to conditional transitivity.
            float tmp = std::min(reference.metric, entry.metric);
            //TransCacheMap::const_iterator peerEntry =
            //    peerInfo->find(entry.estIdx);
            //if (peerEntry == peerInfo->end() ||
            //    peerEntry->second.metric > tmp) {
                (*peerInfo)[entry.estIdx] =
                    CachedESTInfo(estIdx, entry.estIdx,
                                  tmp,
                                  entry.alignmentData);
                //}                
        }
    }
}

bool
TransCacheEntry::getMetric(const int otherESTidx, float& metric) const {
    // First see if our peerInfo hash map has an entry for otherESTidx.
    TransCacheMap::const_iterator peerEntry = peerInfo->find(otherESTidx);
    if (peerEntry == peerInfo->end()) {
        // No, we don't have an entry. Can't determine metric using
        // transitivity in this case.
        return false;
    }
    // OK, found a peer entry.
    metric = peerEntry->second.metric;
    // Return true to indicate we did find an entry.
    return true;
}

#endif
