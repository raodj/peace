#ifndef PRIMES_EST_ANALYZER_CPP
#define PRIMES_EST_ANALYZER_CPP

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

#include "PrimesESTAnalyzer.h"
#include "HeuristicChain.h"

PrimesESTAnalyzer::PrimesESTAnalyzer() : FWAnalyzer("primes") {
    // Initialize command-line args 
    numFeatures   = 4;
    atPrime       = 71;
    cgPrime       = 113;
    distThresh    = -1;
    wordLen       = -1;
}

void
PrimesESTAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    ESTAnalyzer::addCommandLineArguments(argParser);    
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--pri-anal-features", "Number of features for primes heuristic",
         &numFeatures, ArgParser::INTEGER},
        {"--pri-anal-at", "Prime value for A/T in primes heuristic",
         &atPrime, ArgParser::INTEGER},
        {"--pri-anal-cg", "Prime value for C/G in primes heuristic",
         &cgPrime, ArgParser::INTEGER},
        {"--pri-anal-thresh", "Distance threshold override for similarity",
         &distThresh, ArgParser::FLOAT},
        {"--pri-anal-wordLen", "Word length for position-weights",
         &wordLen, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };    
    // Use a arg parser object to conveniently display common options.
    argParser.addValidArguments(ArgsList);
}

bool
PrimesESTAnalyzer::initialize() {
    if (isInitialized()) {
        // Already initialized...
        return true;
    }
    // Let the base class initialize the heuristics in the chain
    if (!ESTAnalyzer::initialize()) {
        // Error occured when initializing.  This is no good.
        return false;
    }
    // Validate command-line arguments.
    if (numFeatures < 1) {
        std::cerr << "primes: Number of features must be >= 1 "
                  << "(use --pri-anal-features option)\n";
        return false;
    }
    // Setup the prime values to use in our base helper class
    setPrimes(atPrime, cgPrime);
    // Generate and cache metrics for local reads
    ASSERT(estList != NULL);
    cacheFeatures(*estList, numFeatures, wordLen);
    // Everying went well
    return true;
}

int
PrimesESTAnalyzer::setReferenceEST(const EST* estS1) {
    ASSERT ( estS1 != NULL );
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(estS1);
    }    
    // If the user has not specified a threshold, compute the average
    // and standard deviation to determine threshold.
    if (distThresh == -1) {
        // User has not specified a threshold. So compute it, based on
        // average distance and variance.  For this we use a helper
        // method in the PrimesHelper base class (that is capable of
        // doing this computation in a distributed manner).
        float average = 0, deviation = 0;
        computeAvgDistance(*estList, estS1->getID(), numFeatures, wordLen,
                           average, deviation);
        // Compute the distance threshold depending on average and
        // deviation
        distThresh = (average > deviation * 1.5 ? (average - deviation) :
                      average);
        std::cerr << "Distance threshold: " << distThresh << " (average: "
                  << average << ", deviation: " << deviation << ")\n";
        ASSERT( distThresh >= 0 );
    }
    
    // Setup the reference sequence features
    refFeatures = extractFeatures(*estS1, numFeatures, wordLen);
    ASSERT( (int) refFeatures.size() == numFeatures );
    // All done for now.
    return 0;
}

float
PrimesESTAnalyzer::getMetric(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    const FloatVec othFeatures = extractFeatures(*otherEST, numFeatures,
                                                 wordLen);
    if (wordLen > 0) {
        // Return the shortest distance between the 2 reads
        return getMinDistance(refFeatures, othFeatures);
    }
    // Non-position-based metric
    const double  distance     = getDistance(refFeatures, othFeatures);
    // return the distance
    return distance;
}

#endif
