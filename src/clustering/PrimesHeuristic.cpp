#ifndef PRIMES_HELPER_CPP
#define PRIMES_HELPER_CPP

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
// Authors:   Dhananjai M. Rao          raodm@miamiOH.edu
//
//---------------------------------------------------------------------

#include <cmath>
#include "EST.h"
#include "PrimesHeuristic.h"
#include "RuntimeContext.h"
#include "MPIStats.h"

PrimesHeuristic::PrimesHeuristic(HeuristicChain* chain, const std::string& name)
    : Heuristic(name, chain) {
    // Initialize command-line args 
    numFeatures   = 4;
    atPrime       = 71;
    cgPrime       = 113;
    distThresh    = -1;
    topN          = -1;
    numSetRef     = 0;
    totSetRefTime = 0;
    topNper     = -1;
}

PrimesHeuristic::~PrimesHeuristic() {
    // Nothing to be done in the destructor
}

void
PrimesHeuristic::addCommandLineArguments(ArgParser& argParser) {
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--pri-heur-features", "Number of features for primes heuristic",
         &numFeatures, ArgParser::INTEGER},
        {"--pri-heur-at", "Prime value for A/T in primes heuristic",
         &atPrime, ArgParser::INTEGER},
        {"--pri-heur-cg", "Prime value for C/G in primes heuristic",
         &cgPrime, ArgParser::INTEGER},
        {"--pri-heur-thresh", "Distance threshold override for similarity",
         &distThresh, ArgParser::LONG},
        {"--pri-heur-topN", "Restrict to top-n reads below threshold",
         &topN, ArgParser::INTEGER},
        {"--pri-heur-topN-per", "Restrict to top-n% of reads below threshold",
         &topNper, ArgParser::FLOAT},
        {"", "", NULL, ArgParser::INVALID}
    };    
    // Use a arg parser object to conveniently display common options.
    argParser.addValidArguments(ArgsList);
}

bool
PrimesHeuristic::initialize() {
    // Validate command-line arguments.
    if (numFeatures < 1) {
        std::cerr << heuristicName << ": Number of features must be >= 1 "
                  << "(use --pri-heur-features option)\n";
        return false;
    }
    // Setup the prime values to use in our base helper class
    setPrimes(atPrime, cgPrime);
    // Generate and cache metrics for local reads
    ESTList& estList = *getContext()->getESTList();
    cacheFeatures(estList, numFeatures);
    // Setup topN if topNper is specified.
    if ((topN == -1) && (topNper > 0)) {
        topN = estList.size() * topNper;
    }
    // Everying went well
    return true;
}

int
PrimesHeuristic::setReferenceEST(const EST* estS1) {
    ASSERT ( estS1 != NULL );
    const double startTime = MPI_Wtime();  // Record the start time
    numSetRef++;  // Track number of times this method is called.
    // Get the EST list to be processed
    ESTList& estList = *getContext()->getESTList();        
    // If the user has not specified a threshold, compute the average
    // and standard deviation to determine threshold.
    if (distThresh == -1) {
        // User has not specified a threshold. So compute it, based on
        // average distance and variance.  For this we use a helper
        // method in the PrimesHelper base class (that is capable of
        // doing this computation in a distributed manner).
        float average = 0, deviation = 0;
        computeAvgDistance(estList, estS1->getID(), numFeatures,
                           average, deviation);
        // Compute the distance threshold depending on average and
        // deviation
        distThresh = (average > deviation * 1.5 ? (average - deviation) :
                      average);
        std::cerr << "Distance threshold: " << distThresh << " (average: "
                  << average << ", deviation: " << deviation << ")\n";
        ASSERT( distThresh > 0 );
    }
    
    // Setup the reference sequence features
    refFeatures = extractFeatures(*estS1, numFeatures);
    ASSERT( (int) refFeatures.size() == numFeatures );

    // Check to see if we need to compute the topN reads if user has
    // specified limits.
    if (topN != -1) {
        // Compute the list of top-n closest ESTs on this node and use
        // that information for further analysis.
        const std::vector<PrimesHelper::ESTMetric> dists =
            computeMetrics(estList, estS1->getID(), numFeatures, distThresh);
        // Restrict to the top-N reads (if set) or dists.size() whichever
        // is smaller.
        const int maxReads = std::min<int>(dists.size(),
                                           (topN != -1) ? topN : dists.size());
        ASSERT(maxReads <= (int) dists.size());
        // Convert the top-n entries into a hash-map for quick look-up
        // in the runHeuristics method.
        nearest.clear();
        for (int i = 0; (i < maxReads); i++) {
            // Add entry to nearest hash map
            nearest[dists[i].estIdx] = dists[i];
        }
    }
    // Track time taken for this call to complete.
    totSetRefTime += (MPI_Wtime() - startTime);
    // All done for now.
    return 0;
}

bool
PrimesHeuristic::runHeuristic(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    // If we have topN set then we can use the nearest hash
    // map. Otherwise we compute the distance and compare against
    // threshold.
    if (topN != -1) {
        // If entry is in the nearest hash map return true. Otherwise
        // return false.
        const int estIdx = otherEST->getID();
        return (nearest.find(estIdx) != nearest.end());
    }
    // In this case we are not using topN but just distance-based cut
    // off
    ASSERT(nearest.empty());
    const LongVec othFeatures = extractFeatures(*otherEST, numFeatures);
    const double  distance    = getDistance(refFeatures, othFeatures);
    // Return true if distance is below threshold
    return (distance <= distThresh);
}

void
PrimesHeuristic::printStats(std::ostream& os) const {
    // First let the base class print common stats
    Heuristic::printStats(os);
    // Print stats from this class
    os << "\tNum. of setRefEST calls: " << numSetRef << "\n"
       << "\tTot time for setRefES  : " << totSetRefTime << " seconds\n";
}

#endif
