#ifndef TWOPASS_D2_CPP
#define TWOPASS_D2_CPP

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

#include "TwoPassD2.h"
#include "EST.h"
#include "ESTCodec.h"
#include "ParameterSetManager.h"
#include <algorithm>

// The bitmask to be used when build hash values.
int TwoPassD2::BitMask   = 0;

// The special value used for words containing an 'N'.
int TwoPassD2::NHash     = 0;

// Instance variable to store the number of bits to be shifted to
// create hash values. This value is initialized to 2*(wordSize-1)
int TwoPassD2::bitShift  = 0;

// The threshold score below which two ESTs are considered
// sufficiently similar to be clustered.
int TwoPassD2::minThreshold = 0;

// Indicates whether or not d2 scores should be normalized based on the
// threshold value.
bool TwoPassD2::noNormalize = false;

// The set of arguments for this class.
arg_parser::arg_record TwoPassD2::argsList[] = {
    {"--d2Threshold", "Threshold score to break out of D2 (default=0)",
     &TwoPassD2::minThreshold, arg_parser::INTEGER},
    {"--noNormalize", "Signals that threshold scores should not be normalized",
     &TwoPassD2::noNormalize, arg_parser::BOOLEAN},
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

TwoPassD2::TwoPassD2(const int refESTidx, const std::string& outputFileName)
    : FWAnalyzer("twopassD2", refESTidx, outputFileName) {
    s1WordTable     = NULL;
    s2WordTable     = NULL;
    delta           = NULL;
    alignmentMetric = 0;
    
    // Below are formerly static parameters that are now dynamic.
    // Note that the following defaults are meaningless as they will be
    // changed before the first comparison is made.
    frameShift   = 50;
    maxThreshold = 130;
    threshold    = 40;
    // Set up current adaptive parameter set index to an invalid value.
    currParamSetIndex = -1;
}

TwoPassD2::~TwoPassD2() {
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
TwoPassD2::showArguments(std::ostream& os) {
    FWAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(TwoPassD2::argsList);
    os << ap;
}

bool
TwoPassD2::parseArguments(int& argc, char **argv) {
    arg_parser ap(TwoPassD2::argsList);
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
TwoPassD2::initialize() {
    // Let the base class initialize any additional heuristics
    int result = FWAnalyzer::initialize();
    if (result != 0) {
        // Error occured when initializing.  This is no good.
        return result;
    }
    // Next let the parameter set manager perform its initialization
    // and optimization of parameter look-ups.
    ParameterSetManager::getParameterSetManager()->initialize();
    // Setup the frequency delta table
    const int MapSize = (1 << (wordSize * 2));
    delta = new int[MapSize];
    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;
    NHash = MapSize;
    // Compute the number of bits to shift when building hashes
    bitShift = 2 * (wordSize - 1);
    // Compute word table size and initialize word tables
    const int wordTableSize = EST::getMaxESTLen() +
        ParameterSetManager::getParameterSetManager()->getMaxFrameSize();
    s1WordTable = new int[wordTableSize];
    s2WordTable = new int[wordTableSize];
    
    // Set the number of words in a window.
    numWordsInWindow = frameSize - wordSize + 1;
    
    // Everything went on well.
    return 0;
}

int
TwoPassD2::setReferenceEST(const int estIdx) {
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(estIdx);
    }
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        // init ref-est word table
        const EST *estS1 = EST::getEST(refESTidx);
        const char* s1   = estS1->getSequence();
        sq1Len = estS1->getSequenceLength();
        // Create the word table using our encoder.
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s1WordTable, s1, encoder);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

float
TwoPassD2::getMetric(const int otherEST) {
    VALIDATE({
            if (otherEST == refESTidx) {
                return 0; // distance to self will be 0
            }
            if ((otherEST < 0) || (otherEST >= EST::getESTCount())) {
                // Invalid EST index!
                return -1;
            }
        });

    // Access information on the comparison sequence
    const EST *estS2 = EST::getEST(otherEST);
    const char* s2   = estS2->getSequence();
    sq2Len           = estS2->getSequenceLength(); 

    // Update the frame size, threshold etc. according to the algorithms
    updateParameters(otherEST);

    // Build the word table for otherEST depending on normal or
    // reverse complement suggestion using hint UVSampleHeuristic.
    int bestMatchIsRC = 0;
    if (chain != NULL) {
        HeuristicChain::getHeuristicChain()->getHint("D2_DoRC", bestMatchIsRC);
    }
    if (bestMatchIsRC) {
        ESTCodec::RevCompEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, s2, encoder);
    } else {
        ESTCodec::NormalEncoder<bitShift, BitMask> encoder;
        buildWordTable(s2WordTable, s2, encoder);
    }

    int s1Index = 0, s2Index = 0, boundDist = 0;
    float normalize = (float) (noNormalize ? 1 : (1.0 / float(threshold)));
    if (frameShift > 1) {
        // OK. Run the asymmetric D2 algorithm
        float distance = (float) runD2Asymmetric(&s1Index, &s2Index);
        if (distance > maxThreshold) {
            return distance * normalize;
        }
        // Set boundDist
        boundDist = frameSize/2;
    } else {
        // If frameshift is 1 we are only running symmetric D2, with no bounds
        // So set the bound distance to be the ends of each sequence
        boundDist = std::max(sq1Len, sq2Len);
    }

    // Now run the bounded symmetric D2 algorithm
    return normalize * (float) runD2Bounded(s1Index-boundDist, 
                                            s1Index+boundDist+frameSize, 
                                            s2Index-boundDist, 
                                            s2Index+boundDist+frameSize);
}

void
TwoPassD2::updateParameters(const int otherEST) {
    const ParameterSetManager* const paramMgr =
        ParameterSetManager::getParameterSetManager();
    const int paramSetIndex = paramMgr->getParameterSet(refESTidx, otherEST);
    if (paramSetIndex == -1) {
        // These two fragments are very different in lengths. Don't
        // bother comparing them at all.
        return;
    }
    if (paramSetIndex == currParamSetIndex) {
        // The parameters that we are currently using are just
        // fine. No need to update them any further. Analyze the two
        // fragments using the current parameter set.
        return;
    }
    // Update and move to a new set of parameters.
    const ParameterSet* const parameterSet = paramMgr->getParameterSet(paramSetIndex);
    // Assign parameters according to parameter set chosen
    frameSize        = parameterSet->frameSize;
    frameShift       = parameterSet->frameShift;
    threshold        = parameterSet->threshold;
    maxThreshold     = parameterSet->maxThreshold;
    numWordsInWindow = frameSize - wordSize + 1;
    // Save the current parameter index for future short circuiting
    currParamSetIndex = paramSetIndex;
}

float
TwoPassD2::runD2Asymmetric(int* s1MinScoreIdx, int* s2MinScoreIdx) {    
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
        if (w1 != NHash) {
            score += (delta[w1] << 1) + 1;
            delta[w1]++;
        }
        // Process i'th word in EST 2
        const int w2 = s2WordTable[sq2Start + i];
        if (w2 != NHash) {
            score -= (delta[w2] << 1) - 1;
            delta[w2]--;
        }    
    }

    // Precompute iteration bounds for the for-loops below to
    // hopefully save on compuation.
    const int LastWindowInSq1  = (sq1End - wordSize + 1) - numWordsInWindow;
    const int LastWindowInSq2  = (sq2End - wordSize + 1) - numWordsInWindow;
    const int FirstWindowInSq2 = sq2Start + numWordsInWindow - 1;
    // Variable to track the minimum d2 distance observed.
    int minScore   = score;
    int s1Win, s2Win, i;
    for(s1Win = sq1Start; (s1Win < LastWindowInSq1 || s1Win == sq1Start);
        s1Win += frameShift*2) {
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
        // Break if the window on s1 cannot be shifted any further right
        if ((s1Win) >= LastWindowInSq1) {
            break;
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
        if (minScore <= minThreshold) {
            break;
        }
        
        // Check every window in EST #2 against current window in EST
        // #1 by sliding EST #2 window to left.
        for(s2Win = sq2End - wordSize; (s2Win > FirstWindowInSq2); s2Win--) {
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
        if (minScore <= minThreshold) {
            break;
        }
    }
    // Now we have the minimum score to report.
    return (float) minScore;
}

float
TwoPassD2::runD2Bounded(int sq1Start, int sq1End, int sq2Start, int sq2End) {
    // Perform sanity checks on bounds (invalid bounds may be passed in)
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
        if (w1 != NHash) {
            score += (delta[w1] << 1) + 1;
            delta[w1]++;
        }
        // Process i'th word in EST 2
        const int w2 = s2WordTable[sq2Start + i];
        if (w2 != NHash) {
            score -= (delta[w2] << 1) - 1;
            delta[w2]--;
        }    
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
    for(s1Win = sq1Start; (s1Win < LastWindowInSq1 || s1Win == sq1Start);
        s1Win += 2) {
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
        // Break if the window on s1 cannot be shifted any further right
        if (s1Win >= LastWindowInSq1) {
            break;
        }
        // Move onto the next window in EST #1.  In this window at
        // (s1Win + numWordsWin) is moving in, while window at s1Win
        // is moving out as we move from left-to-right in EST #1.
        updateWindow(s1WordTable[s1Win], 
                     s1WordTable[s1Win + numWordsInWindow],
                     score, minScore, s1Win+1-s2Win);
        // Break out of this loop if we have found a potential match
        if (minScore <= minThreshold) {
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
        // Break if the window on s1 cannot be shifted any further right
        if ((s1Win+1) >= LastWindowInSq1) {
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
        if (minScore <= minThreshold) {
            break;
        }
    }
    // Now we have the minimum score to report.
    return (float) minScore;
}

bool
TwoPassD2::getAlignmentData(int &alignmentData) {
    // Simply copy the alignment metric that was computed by the last
    // successful call to the analyze() method.
    alignmentData = alignmentMetric;
    // Let the caller know that the alignment data is available.
    return true;
}

#endif
