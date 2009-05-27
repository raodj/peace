#ifndef UV_SAMPLE_HEURISTIC_CPP
#define UV_SAMPLE_HEURISTIC_CPP

//---------------------------------------------------------------------------
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

#include "UVSampleHeuristic.h"
#include "HeuristicChain.h"
#include "ESTCodec.h"
#include "EST.h"

// default params
int UVSampleHeuristic::u           = 4;
int UVSampleHeuristic::v           = 8;
int UVSampleHeuristic::wordShift   = 16;
int UVSampleHeuristic::BitMask     = 0;

// Instance variable to track number of bits to shift for building
// hash values.  This is set to 2*(v-1) in initialize method.
int UVSampleHeuristic::bitsToShift = 0;

// The set of arguments for this class.
arg_parser::arg_record UVSampleHeuristic::argsList[] = {
    {"--uv_u", "u (number of v-word matches) (default=4)",
     &UVSampleHeuristic::u, arg_parser::INTEGER},
    {"--uv_v", "v (length of common words) (default=8)",
     &UVSampleHeuristic::v, arg_parser::INTEGER},
    {"--uv_wordShift", "Word Shift (default=16)",
     &UVSampleHeuristic::wordShift, arg_parser::INTEGER},
    {NULL, NULL}
};

UVSampleHeuristic::UVSampleHeuristic(const std::string& name,
                                     const int refESTIdx,
                                     const std::string& UNREFERENCED_PARAMETER(outputFileName))
    : Heuristic(name, refESTIdx), hintKey("D2_DoRC") {
    // Initialize hash table arrays
    s1WordMap   = NULL;
    s1RCWordMap = NULL;
}

UVSampleHeuristic::~UVSampleHeuristic() {
    if (s1WordMap != NULL) {
        delete [] s1WordMap;
    }
    if (s1RCWordMap != NULL) {
        delete [] s1RCWordMap;
    }
    // Clear out all the entires in the uvCache as they are no longer
    // needed.
    uvCache.clear();
}

void
UVSampleHeuristic::showArguments(std::ostream& os) {
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(UVSampleHeuristic::argsList);
    os << ap;
}

bool
UVSampleHeuristic::parseArguments(int& argc, char **argv) {
    arg_parser ap(UVSampleHeuristic::argsList);
    ap.check_args(argc, argv, false);
    // Ensure values are valid.
    if (wordShift < 1) {
        std::cerr << heuristicName
                  << ": Word shift must be >= 1"
                  << "(use --uv_wordShift option)\n";
        return false;
    } else if (u < 1 || v < 1) {
        std::cerr << heuristicName
                  << ": u and v must be >= 1"
                  << "(use --uv_u and --uv_v options)\n";
        return false;
    }
    // Everything went well
    return true;
}

int
UVSampleHeuristic::initialize() {
    const int MapSize = (1 << (v * 2));
    // Create the normal encoded word map
    s1WordMap = new char[MapSize];
    // Create the reverse-complement encoded word map
    s1RCWordMap = new char[MapSize];
    // Setup the bits to be shifted to create hash values.
    bitsToShift = 2 * (v - 1);
    // Everything went well
    return 0;
}

int
UVSampleHeuristic::setReferenceEST(const int estIdx) {
    if (refESTidx == estIdx) {
        // This EST is already the reference EST. Nothing further to
        // do. So just get out.
        return 0;
    }
    
    // The bit mask has been defined to be a static to enable passing
    // it as a parameter to the templatized NormalEncoder.
    BitMask = (1 << (v * 2)) - 1;
    
    if ((estIdx < 0) || (estIdx >= EST::getESTCount())) {
        // Invalid est index.
        return 1;
    }
    // Setup the look up hash table for the reference est.
    refESTidx = estIdx;
    // Initialize s1 word map for the new reference EST
    const int MapSize = (1 << (v * 2));

    // Initialize word maps to false throughout.
    memset(s1WordMap,   0, sizeof(char) * MapSize);
    memset(s1RCWordMap, 0, sizeof(char) * MapSize);
    // Obtain EST sequence for quick access.
    const EST *estS1 = EST::getEST(refESTidx);
    ASSERT ( estS1 != NULL );
    const char* s1 = estS1->getSequence();
    
    // First compute the hash for a single word using a suitable
    // generator.
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    ESTCodec& codec = ESTCodec::getCodec();
    
    int hash  = 0;
    for(int i = 0; (i < v - 1); i++) {
        hash = encoder(hash, s1[i]);
    }
    // Fill in the word map. But first setup the table to quickly
    // generate reverse complement values.
    codec.setRevCompTable(v);
    const int End = strlen(s1) - v;
    for (int i = v - 1; (i <= End); i++) {
        hash = encoder(hash, s1[i]);
        // Setup the word map for the regular word
        s1WordMap[hash] = 1;
        // Setup the word map entry for the reverse-complement word.
        s1RCWordMap[codec.encode2rc(hash)] = 1;
    }
    return 0; // everything went well
}

void
UVSampleHeuristic::computeHash(const int estIdx) {
    // Pre-compute values that will be used in the for-loop below to
    // reduce duplicate computations.
    const EST *estS2 = EST::getEST(estIdx);
    ASSERT ( estS2  != NULL );
    const char *sq2  = estS2->getSequence();
    ASSERT ( sq2 != NULL );
    const int End    = strlen(sq2) - v;
    ASSERT ( End > 0 );
    
    // Get the codec for encoding/decoding operations

    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    // Obtain the vector from the hash_map and ensure it has
    // sufficient space to hold all the hash values.
    std::vector<unsigned short>& hashList = uvCache[estIdx];
    hashList.reserve(End / wordShift + 1);
    
    // Now, go through the EST sq2 and build hash values
    for (register int start = 0; (start < End); start += wordShift) {
        // Compute the hash for the next v words.
        register unsigned short hash = 0;
        const int endBP    = start + v;
        for(register int i = start; (i < endBP); i++) {
            hash = (unsigned short) encoder(hash, sq2[i]);
        }
        // Store the hash value
        hashList.push_back(hash);
    }
}

bool
UVSampleHeuristic::runHeuristic(const int otherEST) {
    // Extra sanity checks on uncommon scenarios.
    VALIDATE({
        if (otherEST == refESTidx) {
            return true; // will end up with distance 0, or max similarity
        }
        if ((otherEST < 0) || (otherEST >= EST::getESTCount())) {
            // Invalid est index.
            return false;
        }
    });

    // Get otherEST's hash values from the uvCache. If an entry for
    // otherEST is not present in uvCache then build one.
    UVHashTable::iterator cacheEntry = uvCache.find(otherEST);
    if (cacheEntry == uvCache.end()) {
        // An entry does not exist. Create one and cache it for future
        // references.
        computeHash(otherEST);
        // Update our iterator for use below
        cacheEntry = uvCache.find(otherEST);
    }
    // Obtain a reference to the hash list from the cache.
    const std::vector<unsigned short>& otherHash = cacheEntry->second;
    const int hashSize = otherHash.size();
    ASSERT ( hashSize > 0 );
    // go through the otherHash and track number of matching words
    // Initialize local variables.
    register int numMatches = 0, numRCmatches = 0;
    for (register int word = 0; (word < hashSize); word++) {
        // Track number of positive and reverse-complement matches on
        // this hash word
        numMatches   += s1WordMap  [otherHash[word]];
        numRCmatches += s1RCWordMap[otherHash[word]];
    }

    // Print some information for analysis purposes.
    // printf("uv(%d, %d) = %d, %d\n", refESTidx, otherEST, numMatches,
    //       numRCmatches);
    
    // Check if we had sufficient matches to warrant further
    // operations
    if ((numMatches >= u) || (numRCmatches >= u)) {
        // Set the flag to indicate if the normal or the reverse
        // complement version of checks yielded the best result
        bestMatchIsRC = (numMatches < numRCmatches);
        // Setup a hint for D2.
        HeuristicChain::getHeuristicChain()->setHint(hintKey, bestMatchIsRC);
        // return success indication
        return true;
    }
    // When control drops here that means number of matches were not met
    return false;
}

#endif
