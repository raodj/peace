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
#include <numeric>
#include <algorithm>
#include "PrimesHelper.h"
#include "MPIHelper.h"
#include "MPIStats.h"

void
PrimesHelper::setPrimes(const int atPrime, const int cgPrime,
                        const int others) {
    // Ensure we have exactly 256 entries in the numericVector
    numericVector.clear();
    numericVector.resize(256, 0);
    // Set values for specific entries
    numericVector['A'] = numericVector['a'] = +atPrime;
    numericVector['T'] = numericVector['t'] = -atPrime;
    numericVector['U'] = numericVector['u'] = -atPrime;
    numericVector['c'] = numericVector['g'] = +cgPrime;    
    numericVector['C'] = numericVector['c'] = -cgPrime;
    // Convenience list of characters whose values are to be updated
    const std::string Bases = "ryswkmbdhvnRYSWKMBDHVN";
    for (size_t i = 0; (i < Bases.length()); i++) {
        numericVector[i] = others;
    }
}

LongVec
PrimesHelper::extractFeatures(const std::string& sequence,
                              const int numFeatures) const {
    LongVec features(numFeatures);  // The n-dimensional feature
    const int seqLen = sequence.length();  // convenience reference below.
    // Compute scores for each of the n-dimensions.
    for (int startIdx = 0, i = 0; (i < numFeatures); i++) {
        // Convert it to an index position
        const int endIdx = std::min<int>(seqLen, round(seqLen * (i + 1) /
                                                       numFeatures));
        // Compute score.
        features[i] = scoreSequence(sequence, startIdx, endIdx);
        // Update start for next sub-sequence.
        startIdx    = endIdx;
    }
    // Return the computed n-dimensional features back.
    return features;
}

long
PrimesHelper::scoreSequence(const std::string& sequence, const int start,
                            const int end) const {
    long score = 0;
    for (int i = start; (i < end); i++) {
        // The numeric vector for the given character contains the
        // score associated with the given base pair.
        score += numericVector[sequence[i]];
    }
    return score;
}

void
PrimesHelper::cacheFeatures(ESTList& estList, const int features) {
    // First, determine the subset of ESTs to be processed by this
    // MPI-process
    int startIdx, endIdx;
    getLocallyOwnedESTidx(estList.size(), startIdx, endIdx);
    // Compute distances
    for (int idx = startIdx; (idx < endIdx); idx++) {
        const EST* est = estList.get(idx, true);
        featuresCache[idx] = extractFeatures(est->getSequenceString(),
                                             features);
    }
}

// Return value from features cache
LongVec
PrimesHelper::getCachedFeature(const int estIdx) const {
    // Get the iterator corresponding to the estIdx
    FeaturesMap::const_iterator entry = featuresCache.find(estIdx);
    // If we found the features for the given est then return it
    if (entry != featuresCache.end()) {
        return entry->second;  // return pre-computed features
    }
    // Entry not found. Return empty vector
    return LongVec();
}

void
PrimesHelper::computeAvgDistance(ESTList& estList, const int refESTidx,
                                 const int features,
                                 float& average, float& deviation) const {
    // Get the total number of ESTs to be processed.
    const int numESTs = estList.size();
    // First compute the n-dimensional features for the reference EST
    const EST* refEST = estList.get(refESTidx, true);
    const LongVec refFeatures = extractFeatures(refEST->getSequenceString(),
                                                features);
    // Compute the distances between the subset of ESTs owned by this
    // process.  Determine the subset of ESTs to be processed by this
    // MPI-process
    int startIdx, endIdx;
    getLocallyOwnedESTidx(numESTs, startIdx, endIdx);
    // Compute statistics to determine average & variance
    std::vector<double> dists;
    dists.reserve(endIdx - startIdx);
    for (int idx = startIdx; (idx < endIdx); idx++) {
        // See if we have a cached feature. If not compute features
        const LongVec idxFeature = extractFeatures(*estList.get(idx, true),
                                                   features);
        // Compute distnace between reference feature(s)
        const double dist = getDistance(refFeatures, idxFeature);
        // Store into the local dists vector
        dists.push_back(dist);
    }
    // Now compute the sum of the distances (to determine overall average)
    const float localSum = std::accumulate(dists.begin(), dists.end(), 0);
    // Compute global sum of the distances
    float sum = 0;
    MPI_ALL_REDUCE(&localSum, &sum, 1, MPI_TYPE_FLOAT, MPI_OP_SUM);
    average = sum / numESTs;
    // Compute sum of differences to determine variance.
    float localVariance = 0;
    for (size_t i = 0; (i < dists.size()); i++) {
        const float diff = dists[i] - average;
        localVariance += (diff * diff);
    }
    // Compute total sum of differences to determine variance.
    float varSum = 0;
    MPI_ALL_REDUCE(&localVariance, &varSum, 1, MPI_TYPE_FLOAT, MPI_OP_SUM);
    deviation = sqrt(varSum / numESTs);
}

float
PrimesHelper::getDistance(const LongVec& refFeatures,
                          const LongVec& otherFeatures) const {
    // Generate Euclidean distance between the two features.
    ASSERT( otherFeatures.size() == refFeatures.size() );
    const int numFeatures = refFeatures.size();
    long sqDist = 0;
    for (int i = 0; (i < numFeatures); i++) {
        const long diff = refFeatures[i] - otherFeatures[i];
        sqDist += (diff * diff);  // square of distance.
    }
    // The distance is square-root of the sums
    return sqrt(sqDist);
}

LongVec
PrimesHelper::extractFeatures(const EST& est, const int numFeatures) const {
    // Check if entry is in cache.
    LongVec features = getCachedFeature(est.getID());
    // If not compute it.
    if (features.empty()) {
        // Compute the features
        features = extractFeatures(est.getSequenceString(), numFeatures);
    }
    return features;
}

std::vector<PrimesHelper::ESTMetric>
PrimesHelper::computeMetrics(ESTList& estList, int refIdx,
                             int numFeatures) const {
    // Compute reference EST listances
    const LongVec refFeatures = extractFeatures(*estList.get(refIdx, true),
                                                numFeatures);
    // Compute the distances between the subset of ESTs owned by this
    // process.  Determine the subset of ESTs to be processed by this
    // MPI-process
    int startIdx, endIdx;
    getLocallyOwnedESTidx(estList.size(), startIdx, endIdx);
    // Compute statistics to determine average & variance
    std::vector<PrimesHelper::ESTMetric> dists;
    dists.reserve(endIdx - startIdx);
    for (int idx = startIdx; (idx < endIdx); idx++) {
        const EST *est = estList.get(idx, true);
        if (est->hasBeenProcessed()) {
            continue;  // This read can-and-should be ignored.
        }
        // See if we have a cached feature. If not compute features
        const LongVec idxFeature = extractFeatures(*est, numFeatures);
        // Compute distnace between reference feature(s)
        const float dist = getDistance(refFeatures, idxFeature);
        // Store into the local dists vector
        dists.push_back(ESTMetric(idx, dist));
    }
    // Sort the distances and return the distances
    std::sort(dists.begin(), dists.end());
    return dists;
}

#endif
