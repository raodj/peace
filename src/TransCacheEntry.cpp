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
    // Currently we have nothing else to be done here.
}

void
TransCacheEntry::addEntries(const CachedESTInfo& reference,
                            const SMList& metrics,
                            const int UNREFERENCED_PARAMETER(startIndex), 
                            const int UNREFERENCED_PARAMETER(endIndex)) {
    // First ensure that the pre-conditions are met.
    ASSERT (reference.estIdx == estIdx);
    ASSERT (parentInfo.find(reference.refESTidx) == parentInfo.end());
    // Add the reference information to the parent map first.
    parentInfo[reference.refESTidx] = reference;
    // Now add the related metric information to the peerInfo map
    for(size_t index = 0; (index < metrics.size()); index++) {
        const CachedESTInfo& entry = metrics[index];
        ASSERT( entry.refESTidx == reference.refESTidx );
        // Don't add entries referring to ourselves.
        if (entry.estIdx != estIdx) {
            peerInfo[entry.estIdx] = entry;
        }
    }
}

bool
TransCacheEntry::getMetric(const int otherESTidx, float& metric) const {
    // First see if our peerInfo hash map has an entry for
    // otherESTidx.
    TransCacheMap::const_iterator peerEntry = peerInfo.find(otherESTidx);
    if (peerEntry == peerInfo.end()) {
        // No, we don't have an entry. Can't determine metric using
        // transitivity in this case.
        return false;
    }
    // OK, found a peer entry. Obtain the parent entry for this entry
    // from the parentInfo hash map.
    TransCacheMap::const_iterator parentEntry =
        parentInfo.find(peerEntry->second.refESTidx);
    ASSERT ( parentEntry != parentInfo.end() );
    // Now use potential transitivity to compute the metric.
    metric = std::min(parentEntry->second.metric, peerEntry->second.metric);
    // Return true to indicate we did find an entry.
    return true;
}

#endif