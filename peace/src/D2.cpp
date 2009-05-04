#ifndef D2_CPP
#define D2_CPP

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

#include "D2.h"
#include "EST.h"
#include "ESTCodec.h"
#include <algorithm>

// Skip parameter for d2 asymmetric.  Defaults to 1 (i.e. d2 symmetric).
int D2::frameShift = 1;

int D2::BitMask  = 0;
int D2::bitShift = 0;

// The set of arguments for this class.
arg_parser::arg_record D2::argsList[] = {
    {"--frameShift", "Frame Shift (default=1)",
     &D2::frameShift, arg_parser::INTEGER},
    {NULL, NULL}
};

D2::D2(const int refESTidx, const std::string& outputFileName)
    : FWAnalyzer("D2", refESTidx, outputFileName) {
    s1WordTable     = NULL;
    s2WordTable     = NULL;
    delta           = NULL;
    alignmentMetric = 0;
}

D2::~D2() {
    if (s1WordTable != NULL) {
        delete [] s1WordTable;
    }
    if (s2WordTable != NULL) {
        delete [] s2WordTable;
    }
    if (delta != NULL) {
        delete [] delta;
    }
}

void
D2::showArguments(std::ostream& os) {
    FWAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(D2::argsList);
    os << ap;
}

bool
D2::parseArguments(int& argc, char **argv) {
    arg_parser ap(D2::argsList);
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
D2::initialize() {
    // Let the base class initialize any additional heuristics
    int result = FWAnalyzer::initialize();
    if (result != 0) {
        // Error occured when initializing.  This is no good.
        return result;
    }
    // Setup the frequency delta table
    const int MapSize = (1 << (wordSize * 2));
    delta = new int[MapSize];
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;
    // Compute the number of bits to shift when building hashes
    bitShift = 2 * (wordSize - 1);
    // Compute word table size and initialize word tables
    const int wordTableSize = EST::getMaxESTLen() + frameSize;
    s1WordTable = new int[wordTableSize];
    s2WordTable = new int[wordTableSize];
    
    // Set the number of words in a window.
    numWordsInWindow = frameSize - wordSize + 1;
    
    // Everything went on well.
    return 0;
}

int
D2::setReferenceEST(const int estIdx) {
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(estIdx);
    }
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        // init ref-est word table
        const EST *estS1 = EST::getEST(refESTidx);
        const char* s1   = estS1->getSequence();
        // Create the word table using our encoder.
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s1WordTable, s1, encoder);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

float
D2::analyze(const int otherEST) {
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

    // OK. Run the actual d2 algorithm
    return runD2(otherEST);
}

float
D2::runD2(const int otherEST) {
    // Get basic information about the reference EST
    const EST *estS1   = EST::getEST(refESTidx);
    const char* sq1    = estS1->getSequence();
    const int   sq1Len = strlen(sq1);
    // Get basic information about the otherEST EST
    const EST *estS2   = EST::getEST(otherEST);
    const char* sq2    = estS2->getSequence();
    const int   sq2Len = strlen(sq2);
    
    // Build the word table for otherEST depending on normal or
    // reverse complement suggestion using hint UVSampleHeuristic.
    int bestMatchIsRC = 0;
    HeuristicChain::getHeuristicChain()->getHint("D2_DoRC", bestMatchIsRC);
    if (bestMatchIsRC == 0) {
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, sq2, encoder);
    } else {
        ESTCodec::RevCompEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, sq2, encoder);
    }

    // Currently, the bounds on the word compares in d2 is set to the
    // sizes of the two ESTs to compare. However, the bounds can be
    // reduced based on hints from the <i>t/v</i> heuristic.
    int sq1Start = 0, sq1End = sq1Len;
    int sq2Start = 0, sq2End = sq2Len;

    // Initialize the delta tables.
    memset(delta, 0, sizeof(int) * (1 << (wordSize * 2)));
    register int score = 0;
    // First compute the score for first windows.
    for(int i = 0; (i < numWordsInWindow); i++) {
        // Process i'th word in EST 1
        const int w1 = s1WordTable[sq1Start + i];
        score += (delta[w1] << 1) + 1;
        delta[w1]++;
        // Process i'th word in EST 2
        const int w2 = s2WordTable[sq2Start + i];
        score -= (delta[w2] << 1) - 1;
        delta[w2]--;
    }

    // Precompute iteration bounds for the for-loops below to
    // hopefully save on compuation.
    const int LastWordInSq1  = (sq1End - wordSize + 1) - numWordsInWindow;
    const int LastWordInSq2  = (sq2End - wordSize + 1) - numWordsInWindow;
    const int FirstWordInSq2 = sq2Start + numWordsInWindow - 1;
    // Variable to track the minimum d2 distance observed.
    register int minScore   = score;
    for(int s1Win = sq1Start; (s1Win < LastWordInSq1); s1Win += 2) {
        // Check each window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to right
        for(int s2Win = sq2Start; (s2Win < LastWordInSq2); s2Win++) {
            // The word at s2Win + numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from left-to-right.
            updateWindow(s2WordTable[s2Win + numWordsInWindow],
                         s2WordTable[s2Win], score, minScore);
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin) is moving in, while window at s1Win
        // is moving out as we move from left-to-right in EST #1.
        updateWindow(s1WordTable[s1Win], s1WordTable[s1Win + numWordsInWindow],
                     score, minScore);
        // Break out of this loop if we have found a a potential match
        if (minScore <= 40) {
            break;
        }
        
        // Check every window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to left.
        for(int s2Win = sq2End - wordSize; (s2Win > FirstWordInSq2);
            s2Win--) {
            // The word at s2Win - numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from right-to-left.
            updateWindow(s2WordTable[s2Win - numWordsInWindow],
                         s2WordTable[s2Win], score, minScore);
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin + 1) is moving in, while window at
        // s1Win + 1 is moving out as we move from left-to-right in
        // EST #1.
        updateWindow(s1WordTable[s1Win + 1],
                     s1WordTable[s1Win + 1 + numWordsInWindow],
                     score, minScore);
        // Break out of this loop if we have found a a potential match
        if (minScore <= 40) {
            break;
        }
    }
    // Now we have the minimum score to report.
    // printf("d2(%d, %d, %s) = %d\n", refESTidx, otherEST,
    //        bestMatchIsRC ? "rc" : "norm", minScore);
    return minScore;
}

bool
D2::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
