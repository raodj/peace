#ifndef OLD_TWOPASS_D2_CPP
#define OLD_TWOPASS_D2_CPP

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

#include "OldTwoPassD2.h"
#include "ESTList.h"
#include "ESTCodec.h"
#include "ArgParser.h"
#include "HeuristicChain.h"
#include "ParameterSetManager.h"
#include <algorithm>

// The following two are static variables so that they can be passed
// as template parameters (and they need to have external static
// linkage).
// The bitmak to be used when build hash values.
int OldTwoPassD2::BitMask   = 0;
// Instance variable to store the number of bits to be shifted to
// create hash values. This value is initialized to 2*(wordSize-1)
int OldTwoPassD2::bitShift  = 0;

OldTwoPassD2::OldTwoPassD2() : FWAnalyzer("twopassD2") {
    s1WordTable     = NULL;
    s2WordTable     = NULL;
    delta           = NULL;
    // Setup default parameter values.
    frameShift      = 50;
    threshold       = 0;
    maxThreshold    = 130;
    alignmentMetric = 0;
}

OldTwoPassD2::~OldTwoPassD2() {
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
OldTwoPassD2::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    FWAnalyzer::addCommandLineArguments(argParser);
    // Now setup arguments for this class
    const ArgParser::ArgRecord LocalArgsList[] = {
        {"--frameShift", "Frame Shift/skip",
         &frameShift, ArgParser::INTEGER},
        {"--threshold", "Threshold score to break out of D2",
         &threshold, ArgParser::INTEGER},    
        {"--maxThreshold", "Threshold score to run bounded symmetric D2",
         &maxThreshold, ArgParser::INTEGER},    
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add valid arguments to global arg parser
    argParser.addValidArguments(LocalArgsList);
}

bool
OldTwoPassD2::initialize() {
    // Setup parameters so that the heuristics can use them if needed.
    chain->getParamSetMgr()->setupParameters(false);
    // Setup the fact that we would like all 'N' characters to be
    // normalized for this analyzer.
#ifndef _WINDOWS
#warning "Check to ensure randomization of 'N' bases is not yet in place"
    // ESTAnalyzer::setRandomizeNbases(true);
#endif

    // Let the base class initialize any additional heuristics
    if (!FWAnalyzer::initialize()) {
        // Error occured when initializing.  This is no good.
        return false;
    }
    if (frameShift < 1) {
        std::cerr << getName()
                  << ": Frame shift must be >= 1"
                  << "(use --frameShift option)\n";
        return false;
    }
    // initialize the parameter set manager.
    chain->getParamSetMgr()->initialize();
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
    ASSERT( ESTAnalyzer::estList != NULL );
    const int wordTableSize = estList->getMaxESTLen() + frameSize;
    s1WordTable = new int[wordTableSize];
    s2WordTable = new int[wordTableSize];
    
    // Set the number of words in a window.
    numWordsInWindow = frameSize - wordSize + 1;
    
    // Everything went on well.
    return true;
}

int
OldTwoPassD2::setReferenceEST(const EST* est) {
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(est);
    }
    // Update the reference EST to be used
    refEST = est;
    if (est != NULL) {
        // init ref-est word table
        const char* s1   = refEST->getSequence();
        // Create the word table using our encoder.
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s1WordTable, s1, encoder);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

float
OldTwoPassD2::getMetric(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    VALIDATE({
            if (otherEST->getID() == refESTidx->getID()) {
                return 0; // distance to self will be 0
            }
        });

    int s1Index = 0;
    int s2Index = 0;
    
    // OK. Run the asymmetric D2 algorithm
    float distance = (float) runD2Asymmetric(otherEST, &s1Index, &s2Index);
    if (distance > maxThreshold) {
        return distance;
    }
    else {
        // Now run the bounded symmetric D2 algorithm
        int boundDist = frameShift/2;
        return (float) runD2Bounded(otherEST, s1Index-boundDist, 
                                    s1Index+boundDist+frameSize, 
                                    s2Index-boundDist, 
                                    s2Index+boundDist+frameSize);
    }
}

float
OldTwoPassD2::runD2Asymmetric(const EST* otherEST, int* s1MinScoreIdx, 
                              int* s2MinScoreIdx) {
    // Get basic information about the reference EST
    const int   sq1Len = refEST->getSequenceLength();
    // Get basic information about the otherEST
    const char* sq2    = otherEST->getSequence();
    const int   sq2Len = otherEST->getSequenceLength();
    
    // Build the word table for otherEST depending on normal or
    // reverse complement suggestion using hint UVSampleHeuristic.
    int bestMatchIsRC = 0;
    if (chain != NULL) {
        chain->getHint(HeuristicChain::D2_DO_RC, bestMatchIsRC);
    }
    if (bestMatchIsRC) {
        ESTCodec::RevCompEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, sq2, encoder);
    } else {
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, sq2, encoder);
    }

    // Currently, the bounds on the word compares in d2 is set to the
    // sizes of the two ESTs to compare. However, the bounds can be
    // reduced based on hints from the <i>t/v</i> heuristic.
    int sq1Start = 0, sq1End = sq1Len;
    int sq2Start = 0, sq2End = sq2Len;

    // Initialize the delta tables.
    memset(delta, 0, sizeof(int) * (1 << (wordSize * 2)));
    int score = 0;
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
    const int LastWindowInSq1  = (sq1End - wordSize + 1) - numWordsInWindow;
    const int LastWindowInSq2  = (sq2End - wordSize + 1) - numWordsInWindow;
    const int FirstWindowInSq2 = sq2Start + numWordsInWindow - 1;
    // Variable to track the minimum d2 distance observed.
    int minScore   = score;
    int s1Win, s2Win, i;
    for(s1Win = sq1Start; (s1Win < LastWindowInSq1); s1Win += frameShift*2) {
        // Check each window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to right
        for(s2Win = sq2Start; (s2Win < LastWindowInSq2); s2Win++) {
            // The word at s2Win + numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from left-to-right.
            updateWindowAsym(s2WordTable[s2Win + numWordsInWindow],
                             s2WordTable[s2Win], score, minScore, 
                             s1Win, s2Win+1, s1MinScoreIdx, s2MinScoreIdx);
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin) is moving in, while window at s1Win
        // is moving out as we move from left-to-right in EST #1.
        // For asymmetric D2, we shift s1 more than one word at a time
        for (i = 0; (i < frameShift && s1Win+i < LastWindowInSq1); i++) {
            updateWindowAsym(s1WordTable[s1Win + i], 
                             s1WordTable[s1Win + i + numWordsInWindow],
                             score, minScore, 
                             s1Win+i+1, LastWindowInSq2,
                             s1MinScoreIdx, s2MinScoreIdx);
        }
        // Break out of this loop if we have found a potential match
        if (minScore <= threshold) {
            break;
        }
        
        // Check every window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to left.
        for(s2Win = sq2End - wordSize; (s2Win > FirstWindowInSq2);
            s2Win--) {
            // The word at s2Win - numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from right-to-left.
            updateWindowAsym(s2WordTable[s2Win - numWordsInWindow],
                             s2WordTable[s2Win], score, minScore, 
                             s1Win+i, s2Win-numWordsInWindow, 
                             s1MinScoreIdx, s2MinScoreIdx);
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin + 1) is moving in, while window at
        // s1Win + 1 is moving out as we move from left-to-right in
        // EST #1.
        // For asymmetric D2, we shift s1 more than one word at a time
        for (i = frameShift; (i < frameShift*2 && s1Win+i < LastWindowInSq1);
             i++) {
            updateWindowAsym(s1WordTable[s1Win + i],
                             s1WordTable[s1Win + i + numWordsInWindow],
                             score, minScore, 
                             s1Win+i+1, 0, 
                             s1MinScoreIdx, s2MinScoreIdx);
        }
        // Break out of this loop if we have found a potential match
        if (minScore <= threshold) {
            break;
        }
    }
    // Now we have the minimum score to report.
    return (float) minScore;
}

float
OldTwoPassD2::runD2Bounded(const EST* otherEST, int sq1Start, int sq1End, 
                           int sq2Start, int sq2End) {
    // Get basic information about the reference EST
    const int   sq1Len = refEST->getSequenceLength();
    // Get basic information about the otherEST EST
    const int   sq2Len = otherEST->getSequenceLength();
    
    // sanity checks on bounds (invalid bounds may be passed in)
    if (sq1Start < 0) sq1Start = 0;
    if (sq2Start < 0) sq2Start = 0;
    if (sq1End > sq1Len) sq1End = sq1Len;
    if (sq2End > sq2Len) sq2End = sq2Len;
	
    // Initialize the delta tables.
    memset(delta, 0, sizeof(int) * (1 << (wordSize * 2)));
    int score = 0;
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
    const int LastWindowInSq1  = (sq1End - wordSize + 1) - numWordsInWindow;
    const int LastWindowInSq2  = (sq2End - wordSize + 1) - numWordsInWindow;
    const int FirstWindowInSq2 = sq2Start + numWordsInWindow - 1;

    // Variable to track the minimum d2 distance observed.
    int minScore   = score;
    alignmentMetric = sq1Start - sq2Start;
    int s1Win, s2Win;
    for(s1Win = sq1Start; (s1Win < LastWindowInSq1); s1Win += 2) {
        // Check each window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to right
        for(s2Win = sq2Start; (s2Win < LastWindowInSq2); s2Win++) {
            // The word at s2Win + numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from left-to-right.
            updateWindow(s2WordTable[s2Win + numWordsInWindow],
                         s2WordTable[s2Win], score, minScore,
                         s1Win-s2Win-1);
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin) is moving in, while window at s1Win
        // is moving out as we move from left-to-right in EST #1.
        updateWindow(s1WordTable[s1Win], 
                     s1WordTable[s1Win + numWordsInWindow],
                     score, minScore, s1Win+1-s2Win);
        // Break out of this loop if we have found a potential match
        if (minScore <= threshold) {
            break;
        }
        
        // Check every window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to left.
        for(s2Win = sq2End - wordSize; (s2Win > FirstWindowInSq2); s2Win--) {
            // The word at s2Win - numWordsInWindow is moving in while
            // the word at s2Win is moving out as we move window
            // associated with EST #2 from right-to-left.
            updateWindow(s2WordTable[s2Win - numWordsInWindow],
                         s2WordTable[s2Win], score, minScore,
                         s1Win+1-s2Win+numWordsInWindow);
        }
        if ((s1Win+1) == LastWindowInSq1) {
            break;
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin + 1) is moving in, while window at
        // s1Win + 1 is moving out as we move from left-to-right in
        // EST #1.
        updateWindow(s1WordTable[s1Win + 1],
                     s1WordTable[s1Win + 1 + numWordsInWindow],
                     score, minScore, s1Win+2-sq2Start);
        // Break out of this loop if we have found a potential match
        if (minScore <= threshold) {
            break;
        }
    }
    // Now we have the minimum score to report.
    return (float) minScore;
}

bool
OldTwoPassD2::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
