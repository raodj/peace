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
#include "EST.h"
#include "ESTCodec.h"

#include <cmath>

// default params
int UVSampleHeuristic::u = 2;
int UVSampleHeuristic::v = 10;
int UVSampleHeuristic::wordShift = 16;
int UVSampleHeuristic::BitMask = 0;


// The set of arguments for this class.
arg_parser::arg_record UVSampleHeuristic::argsList[] = {
    {"--uv_u", "u (number of v-word matches) (default=2)",
     &UVSampleHeuristic::u, arg_parser::INTEGER},
    {"--uv_v", "v (length of common words) (default=10)",
     &UVSampleHeuristic::v, arg_parser::INTEGER},
    {"--uv_wordShift", "Word Shift (default=16)",
     &UVSampleHeuristic::wordShift, arg_parser::INTEGER},
    {NULL, NULL}
};

UVSampleHeuristic::UVSampleHeuristic(const int refESTIdx,
                                     const std::string& outputFileName)
    : Heuristic("uv", refESTIdx) {
}

UVSampleHeuristic::~UVSampleHeuristic() {
    if (s1WordMap != NULL) {
        delete [] s1WordMap;
    }
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
    const int MapSize = pow(4, v);
    s1WordMap = new char[MapSize];
    // Everything went well
    return 0;
}

int
UVSampleHeuristic::setReferenceEST(const int estIdx) {
    // The bit mask has been defined to be a static to enable passing
    // it as a parameter to the templatized NormalEncoder.
    BitMask = (1 << (v * 2)) - 1;
    
    // Codec for encoding/decoding operations
    const ESTCodec& codec = ESTCodec::getCodec();
    
    if ((estIdx < 0) || (estIdx >= EST::getESTCount())) {
        // Invalid est index.
        return 1;
    }
    // Setup the look up hash table for the reference est.
    refESTidx = estIdx;
    // Initialize s1 word map for the new reference EST
    const int MapSize = pow(4, v);

    // Initialize word map to false throughout.
    memset(s1WordMap, 0, sizeof(char) * MapSize);
    // Obtain EST sequence for quick access.
    const EST *estS1 = EST::getEST(refESTidx);
    ASSERT ( estS1 != NULL );
    const char* s1 = estS1->getSequence();
    
    // First compute the hash for a single word using a suitable
    // generator.
    ESTCodec::NormalEncoder<v, BitMask> encoder;
    int hash = 0;
    for(int i = 0; (i < v); i++) {
        hash = encoder(hash, s1[i]);
    }
    // Setup the word map entry for the first word.
    s1WordMap[hash] = 1;
    // Fill in the word map
    const int End = strlen(s1) - v;
    for (int i = 1; (i <= End); i++) {
        hash = encoder(hash, s1[i]);
        s1WordMap[hash] = 1;
    }
    return 0; // everything went well
}

bool
UVSampleHeuristic::runHeuristic(const int otherEST) {
    if (otherEST == refESTidx) {
        return true; // will end up with distance 0, or max similarity
    }
    if ((otherEST < 0) || (otherEST >= EST::getESTCount())) {
        // Invalid est index.
        return false;
    }
    int numMatches = 0;
    // Get s2 sequence (we're done with s1 at this point)
    EST *estS2 = EST::getEST(otherEST);
    std::string sq2 = estS2->getSequence();
    ASSERT ( sq2.size() > 0);

    // Get the codec for encoding/decoding operations
    const ESTCodec& codec = ESTCodec::getCodec();
    const int BitMask     = (1 << (v * 2)) - 1;
    int hash = 0;
    // get initial v-word
    for(int i = 0; (i < v); i++) {
        hash <<= 2;
        hash  |= codec.encode(sq2[i]);
    }
    
    // Track if this word is found in s1
    numMatches += s1WordMap[hash];
    
    // go through the rest of s2 and check
    const int End = sq2.size() - v + 1;
    for (int i = 1; ((i <= End) && (numMatches < u)); i++) {
        hash <<= 2;
        hash  |= codec.encode(sq2[i]);
        hash  &= BitMask;
        numMatches += s1WordMap[hash];
    }
    
    //printf("UV heuristic: %d %d %d\n", refESTidx, otherEST, numMatches);
    return (numMatches >= u);
}

#endif
