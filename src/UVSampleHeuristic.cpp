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
#include <cmath>

#define MAX_SEQ_LEN 5000

char UVSampleHeuristic::CharToInt[256];

// default params
int UVSampleHeuristic::u = 2;
int UVSampleHeuristic::v = 10;
int UVSampleHeuristic::wordShift = 16;

static int BitMask;
static int MapSize;

// Indicates whether or not the given v-word appears in s1
bool* s1WordMap;

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
    // Initialize the array CharToInt for mapping A, T, C, and G to 0,
    // 1, 2, and 3 respectively.
    memset(CharToInt, 0, sizeof(CharToInt));
    // Initialize values for specific base pairs.
    CharToInt[(int) 'A'] = CharToInt[(int) 'a'] = 0;
    CharToInt[(int) 'G'] = CharToInt[(int) 'g'] = 1;
    CharToInt[(int) 'C'] = CharToInt[(int) 'c'] = 2;
    CharToInt[(int) 'T'] = CharToInt[(int) 't'] = 3;
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
    MapSize = pow(4, v);

    BitMask = (1 << (v * 2)) -1;

    s1WordMap = new bool[MapSize];

    // Everything went well
    return 0;
}

int
UVSampleHeuristic::setReferenceEST(const int estIdx) {
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        // init s1 word map
        memset(s1WordMap, 0, sizeof(bool) * MapSize);
        EST *estS1 = EST::getEST(refESTidx);
        std::string s1 = estS1->getSequence();
        // First compute the hash for a single word.
        //ASSERT ( sequence != NULL );
        int hash = 0;
        int i;
        for(i = 0; (i < v); i++) {
            hash <<= 2;
            hash  |= CharToInt[(int) s1[i]];
        }
        // Setup the word map entry for the first word.
        s1WordMap[hash] = true;
        // Fill in the word map
        for (i = 1; (i+v <= s1.size()); i++) {
            hash <<= 2;
            hash  |= CharToInt[(int) s1[i]];
            hash  &= BitMask;
            s1WordMap[hash] = true;
        }
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

bool
UVSampleHeuristic::runHeuristic(const int otherEST) {
    if (otherEST == refESTidx) {
        return true; // will end up with distance 0, or max similarity
    } else if ((otherEST >= 0) && (otherEST < EST::getESTCount())) {
        int numMatches = 0;
        // Get s2 sequence (we're done with s1 at this point)
        EST *estS2 = EST::getEST(otherEST);
        std::string sq2 = estS2->getSequence();
        ASSERT ( sq2.size() > 0);

        // main logic
        int hash = 0;
        int i;
        
        // get initial v-word
        for(i = 0; (i < v); i++) {
            hash <<= 2;
            hash  |= CharToInt[(int) sq2[i]];
        }

        // check if this word is found in s1
        if (s1WordMap[hash]) numMatches++;

        // go through the rest of s2 and check
        for (i = 1; (i+v <= sq2.size()); i++) {
            if (numMatches >= u) break; // have enough matches
            hash <<= 2;
            hash  |= CharToInt[(int) sq2[i]];
            hash  &= BitMask;
            if (s1WordMap[hash]) numMatches++;
        }

        //printf("UV heuristic: %d %d %d\n", refESTidx, otherEST, numMatches);

        return (numMatches >= u);
    } else {
        // Invalid est index
        return false;
    }
}


#endif
