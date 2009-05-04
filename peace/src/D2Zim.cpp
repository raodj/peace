#ifndef D2_ZIM_CPP
#define D2_ZIM_CPP

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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "D2Zim.h"
#include "EST.h"
#include "ESTCodec.h"
#include <algorithm>

// Skip parameter for d2 asymmetric.  Defaults to 1 (i.e. d2 symmetric).
int D2Zim::frameShift = 1;

int D2Zim::BitMask = 0;

// Size of the word tables
size_t D2Zim::wordTableSize = 0;

// The set of arguments for this class.
arg_parser::arg_record D2Zim::argsList[] = {
    {"--frameShift", "Frame Shift (default=1)",
     &D2Zim::frameShift, arg_parser::INTEGER},
    {NULL, NULL}
};

D2Zim::D2Zim(const int refESTidx, const std::string& outputFileName)
    : FWAnalyzer("D2", refESTidx, outputFileName) {
    fdHashMap   = NULL;
    s1WordTable = NULL;
    s2WordTable = NULL;
}

D2Zim::~D2Zim() {
    if (fdHashMap != NULL) {
        delete [] fdHashMap;
    }
    if (s1WordTable != NULL) {
        delete [] s1WordTable;
    }
    if (s2WordTable != NULL) {
        delete [] s2WordTable;
    }
}

void
D2Zim::showArguments(std::ostream& os) {
    FWAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(D2Zim::argsList);
    os << ap;
}

bool
D2Zim::parseArguments(int& argc, char **argv) {
    arg_parser ap(D2Zim::argsList);
    ap.check_args(argc, argv, false);
    // Ensure frameshift is valid.
    if (frameShift < 1) {
        std::cerr << analyzerName
                  << ": Frame shift must be >= 1"
                  << "(use --frameShift option)\n";
        return false;
    }
    // Now let the base class do processing and return the result.
    return FWAnalyzer::parseArguments(argc, argv);
}

int
D2Zim::initialize() {
    int result = FWAnalyzer::initialize();
    if (result != 0) {
        // Error occured when loading ESTs.  This is no good.
        return result;
    }
    // Compute the size of the frequency differential hashmap
    const int MapSize = 1 << (wordSize * 2);

    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;

    NumWordsWin = frameSize-wordSize;
    
    // Initialize the fd hash map

    fdHashMap = new int[MapSize];

    // Compute word table size and initialize word tables
    wordTableSize = EST::getMaxESTLen() + frameSize;
    s1WordTable   = new int[wordTableSize];
    s2WordTable   = new int[wordTableSize];
    
    // Everything went on well.
    return 0;
}

int
D2Zim::setReferenceEST(const int estIdx) {
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(estIdx);
    }
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        // init ref-est word table
        EST *estS1 = EST::getEST(refESTidx);
        const char* s1 = estS1->getSequence();
        buildWordTable(s1WordTable, s1);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

void
D2Zim::buildWordTable(int* wordTable, const char* s) {
    // init ref-est word table
    memset(wordTable, 0, sizeof(int) * wordTableSize);
    // First compute the hash for a single word.
    ASSERT ( s != NULL );
    const int length = strlen(s) - wordSize;
    ESTCodec::NormalEncoder<wordSize, BitMask> encoder;
    register int hash = 0;
    for(int i = 0; (i < wordSize); i++) {
        hash = encoder(hash, s[i]);
    }
    // Setup the word table entry for the first word.
    wordTable[0] = hash;
    // Fill in the word table
    for (int i = 1; (i <= length); i++) {
        hash = encoder(hash, s[i]);
        wordTable[i] = hash;
    }
}

void
D2Zim::buildFdHashMaps(int* sed) {
    // Clear out any old entries in the hash maps.
    const int MapSize = 1 << (wordSize * 2);
    memset(fdHashMap, 0, sizeof(int) * MapSize);

    for (int i = 0; (i+wordSize <= frameSize); i++) {
        *sed +=  2*(fdHashMap[s1WordTable[i]]++) + 1;
        *sed += -2*(fdHashMap[s2WordTable[i]]--) + 1;
    }
}

// This common code fragement is used 4 times in the D2Zim::analyze()
// method below. It was pulled out into a macro to streamline the
// analyze method.
#define CHECK_SED_AND_BREAK           \
    if (sed < minSed) {               \
        alignmentMetric = s1FramePos; \
        if ((minSed = sed) == 0) {    \
            break;                    \
        }                             \
}

float
D2Zim::analyze(const int otherEST) {
    if (otherEST == refESTidx) {
        return 0; // distance to self will be 0
    }
    if ((otherEST < 0) || (otherEST >= EST::getESTCount())) {
        // Invalid EST index!
        return -1;
    }
    // Check with the heuristic chain
    if ((chain != NULL) && (!chain->shouldAnalyze(otherEST))) {
        // Heuristics indicate we should not do D2. So skip it.
        return (4 * frameSize);
    }
    
    // Perform operations for D2
    int sed = 0;
    int minSed = frameSize * 4; // won't be exceeded
    int s1FramePos = 0;
    // Get sequences
    EST *estS1 = EST::getEST(refESTidx);
    EST *estS2 = EST::getEST(otherEST);
    const char* sq1 = estS1->getSequence();
    ASSERT ( sq1 != NULL );
    const char* sq2 = estS2->getSequence();
    ASSERT ( sq2 != NULL );

    // Initialize the word table for s2
    buildWordTable(s2WordTable, sq2);
    
    // Initialize the frequency differential hashmap for s1-s2
    buildFdHashMaps(&sed);
    
    minSed = sed;

    const int numFramesS1 = strlen(sq1) - frameSize;
    const int numFramesS2 = strlen(sq2) - frameSize;
    
    // Main d2 algorithm (from Zimmermann paper)
    while (s1FramePos < numFramesS1) {
        if (s1FramePos != 0) {
            for(int i = 0; (i < frameShift && s1FramePos++ < numFramesS1); i++){
                refShiftUpdateFd(&sed, s1FramePos);
                CHECK_SED_AND_BREAK;
            }
        }
        for(int s2FramePos = 1; s2FramePos <= numFramesS2; s2FramePos++) {
            rightShiftUpdateFd(&sed, s2FramePos);
            CHECK_SED_AND_BREAK;
        }
        if (s1FramePos != numFramesS1) {
            for(int i = 0; (i < frameShift && s1FramePos++ < numFramesS1); i++){
                refShiftUpdateFd(&sed, s1FramePos);
                CHECK_SED_AND_BREAK;
            }
        }
        for(int s2FramePos = numFramesS2-1; s2FramePos >= 0; s2FramePos--) {
            leftShiftUpdateFd(&sed, s2FramePos);
            CHECK_SED_AND_BREAK;
        }
    }
    
    //printf("%d %d %d\n", refESTidx, otherEST, minSed);   
    return minSed;
}

bool
D2Zim::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
