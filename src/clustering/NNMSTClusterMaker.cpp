#ifndef NN_MST_CLUSTER_MAKER_CPP
#define NN_MST_CLUSTER_MAKER_CPP

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

#include "ESTList.h"
#include "MPIStats.h"
#include "HeuristicChain.h"
#include "NNMSTClusterMaker.h"

#include <iostream>

NNMSTClusterMaker::NNMSTClusterMaker(ESTAnalyzer *analyzer)
    : MSTClusterMaker("nnMST", analyzer) {
    // Initialize other instance variables & command-line arguments to
    // custom default values.
    localESTStartIndex = -1;
}

NNMSTClusterMaker::~NNMSTClusterMaker() {
    // Nothing to be done for now.
}

bool
NNMSTClusterMaker::initialize() {
    if (!MSTClusterMaker::initialize()) {
        // Base class inititalization failed
        return false;
    }
    // Get the locally-owned contiguous range of ESTs to setup the
    // nearest-neighbor pool.
    int lastESTindex = -1;
    getLocallyOwnedESTidx(localESTStartIndex, lastESTindex);
    // Setup the nearest-neighbor pool with all false entries.
    nnPool.resize(lastESTindex - localESTStartIndex + 1, false);
    // Everything went well.
    return true;
}


bool
NNMSTClusterMaker::computeSMEntry(const int estIdx,
                                  const int otherEstIdx, SMList& smList,
                                  const float InvalidMetricIgn,
                                  const bool addInvalidMetric,
                                  const bool addToNNPool) {
    // Get similarity/distance metric.
    const float metric = analyze(estList->get(otherEstIdx));
    ASSERT ( metric >= 0 );
    // Obtain any alignment data that the analyzer may have.
    int alignmentData = 0;
    int directionData = 1;
    analyzer->getAlignmentData(alignmentData);
    const HeuristicChain* const chain = analyzer->getHeuristicChain();
    if (chain != NULL) {
        chain->getHint(HeuristicChain::MST_RC, directionData);
    }
    // Add only the first invalid entry. One is enough to do build
    // a valid MST. There is no need for multiple vestigial
    // entries.
    const float InvalidMetric = analyzer->getInvalidMetric();
    const bool isMetricOK = analyzer->compareMetrics(metric, InvalidMetric);
    if ((isMetricOK) || (addInvalidMetric)) {
        // Add the information to metric list
        smList.push_back(CachedESTInfo(estIdx, otherEstIdx,
                                       metric, alignmentData,
                                       directionData));
        // Add valid entries to the nearest-neighbor pool (as directed)
        if (isMetricOK && addToNNPool) {
            nnPool[otherEstIdx - localESTStartIndex] = true;
        }
        // Return true to indicate entry was added.
        return true;
    }
    // Return false to indicate entry was not added
    return false;
}

void
NNMSTClusterMaker::computeSMList(const int estIdx, SMList& smList) {
    // Mark the est that has been added to MST as having been processed
    estList->get(estIdx)->setProcessed(true);
    // Pre-allocate sufficient space in outgoing vector to minimize
    // repeated memory growth. We intentionally do not use endESTidx -
    // startESTidx for size because the number of ESTs can be in
    // millions and will cause memory issues.
    smList.reserve(128);

    printNNPool(std::cout);
    
    // Determine the list of ESTs that this process must deal
    // with using the helper method in base class.
    int startESTidx, endESTidx;
    getLocallyOwnedESTidx(startESTidx, endESTidx);

    // Setup the reference estIdx in the analyzer which given the
    // analyzer a chance to optimize initialization.
    ASSERT( estList != NULL );
    analyzer->setReferenceEST(estList->get(estIdx));
    // Short-cut to invalid metric returned by analyzer based on
    // length of the reference EST.
    const float InvalidMetric = analyzer->getInvalidMetric();

    // First search in the nearest-neighbor pool for a match.
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if (estList->get(otherIdx)->hasBeenProcessed() ||
            !nnPool[otherIdx - startESTidx]) {
            // This EST entry can be ignored because:
            // 1. This metric is not needed (as node is in MST) OR
            // 2. It must not be used (as it has been filtered out) OR
            // 3. It is not in the nearest-neighbor pool.
            continue;
        }
        // Add only valid/successful analysis results in this phase
        computeSMEntry(estIdx, otherIdx, smList, InvalidMetric, false, false);
    }

    // Distribute & obtain the total number of matches found by all
    // the other processes to determine if further operations are
    // needed.
    int localCount = smList.size();
    int totalCount = 0;
    MPI_ALL_REDUCE(&localCount, &totalCount, 1, MPI_TYPE_INT, MPI_OP_SUM);
    // If we found some matches from the nn-pool we don't need to do
    // the next phase.
    if (totalCount > 0) {
        // Found one-or-more entries from the nnPool. Good -- no
        // further searching is needed.
        return;
    }

    // When control drops here, we did not find any matches in the
    // nn-pool. So we need to exhaustively search all locally-owned
    // ESTs.
    
    // Flag to ensure only one invalid metric gets added
    bool needInvalidMetric = true;
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if ((estList->get(otherIdx)->hasBeenProcessed()) ||
            (nnPool[otherIdx - startESTidx])) {
            // This EST entry can be ignored as this metric is not
            // needed (or must not be used) or it has already been
            // checked in the previous phase
            continue;
        }
        // Add similarity/distance metric (as needed)
        if (computeSMEntry(estIdx, otherIdx, smList,
                           InvalidMetric, needInvalidMetric, true)) {
            needInvalidMetric = false;
        }
    }
}

void
NNMSTClusterMaker::printNNPool(std::ostream& os) {
    // Determine the list of ESTs that this nnPool deals with.
    int startESTidx, endESTidx;
    getLocallyOwnedESTidx(startESTidx, endESTidx);

    // Print entries present in the nnPool
    for(int otherIdx = startESTidx; (otherIdx < endESTidx); otherIdx++) {
        if (nnPool[otherIdx - startESTidx]) {
            os << otherIdx << " ";
        }
    }
}


#endif
