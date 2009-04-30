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
#include <cmath>
#include <algorithm>

// The simple translation table to convert chars to int.
char D2::CharToInt[256];

// Skip parameter for d2 asymmetric.  Defaults to 1 (i.e. d2 symmetric).
int D2::frameShift = 1;

static int MapSize;

static int NumWordsWin;

int D2::BitMask = 0;

// Size of the word tables
size_t D2::wordTableSize = 0;

// The set of arguments for this class.
arg_parser::arg_record D2::argsList[] = {
    {"--frameShift", "Frame Shift (default=1)",
     &D2::frameShift, arg_parser::INTEGER},
    {NULL, NULL}
};

D2::D2(const int refESTidx, const std::string& outputFileName)
    : FWAnalyzer("D2", refESTidx, outputFileName) {
    // Initialize the array CharToInt for mapping A, T, C, and G to 0,
    // 1, 2, and 3 respectively.
    memset(CharToInt, 0, sizeof(CharToInt));
    // Initialize values for specific base pairs.
    CharToInt[(int) 'A'] = CharToInt[(int) 'a'] = 0;
    CharToInt[(int) 'G'] = CharToInt[(int) 'g'] = 1;
    CharToInt[(int) 'C'] = CharToInt[(int) 'c'] = 2;
    CharToInt[(int) 'T'] = CharToInt[(int) 't'] = 3;
    fdHashMap = NULL;
    s1WordTable = NULL;
    s2WordTable = NULL;
}

D2::~D2() {
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
    int result = FWAnalyzer::initialize();
    if (result != 0) {
        // Error occured when loading ESTs.  This is no good.
        return result;
    }
    // Compute the size of the frequency differential hashmap
    MapSize = pow(4, wordSize);

    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;

    NumWordsWin = frameSize-wordSize;
    
    // Initialize the fd hash map

    fdHashMap = new int[MapSize];

    // Compute word table size and initialize word tables
    wordTableSize = EST::getMaxESTLen() + frameSize;
    s1WordTable = new int[wordTableSize];
    s2WordTable = new int[wordTableSize];
    
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
        EST *estS1 = EST::getEST(refESTidx);
        std::string s1 = estS1->getSequence();
        buildWordTable(s1WordTable, s1);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

void
D2::buildWordTable(int* wordTable, std::string s) {
    // init ref-est word table
    memset(wordTable, 0, sizeof(int) * wordTableSize);
    // First compute the hash for a single word.
    //ASSERT ( sequence != NULL );
    int hash = 0;
    int i;
    for(i = 0; (i < wordSize); i++) {
        hash <<= 2;
        hash  |= CharToInt[(int) s[i]];
    }
    // Setup the word table entry for the first word.
    wordTable[0] = hash;
    // Fill in the word table
    for (i = 1; (i+wordSize <= (int) s.size()); i++) {
        hash <<= 2;
        hash  |= CharToInt[(int) s[i]];
        hash  &= BitMask;
        wordTable[i] = hash;
    }
}

void
D2::buildFdHashMaps(int* sed) {
    // Clear out any old entries in the hash maps.
    memset(fdHashMap, 0, sizeof(int) * MapSize);

    for (int i = 0; (i+wordSize <= frameSize); i++) {
        *sed +=  2*(fdHashMap[s1WordTable[i]]++) + 1;
        *sed += -2*(fdHashMap[s2WordTable[i]]--) + 1;
    }
}

void
D2::refShiftUpdateFd(int* sed, const int framePos) {
    // update sed and fd from leftmost word falling out
    *sed += -2*(fdHashMap[s1WordTable[framePos-1]]--) + 1;
    // update sed and fd from new rightmost word
    *sed +=  2*(fdHashMap[s1WordTable[framePos+NumWordsWin]]++) + 1;
}

void
D2::rightShiftUpdateFd(int* sed, const int framePos) {
    // update sed and fd from leftmost word falling out
    *sed +=  2*(fdHashMap[s2WordTable[framePos-1]]++) + 1;
    // update sed and fd from new rightmost word
    *sed += -2*(fdHashMap[s2WordTable[framePos+NumWordsWin]]--) + 1;
}

void
D2::leftShiftUpdateFd(int* sed, const int framePos) {
    // update sed and fd from rightmost word falling out
    *sed +=  2*(fdHashMap[s2WordTable[framePos+NumWordsWin+1]]++) + 1;
    // update sed and fd from new leftmost word
    *sed += -2*(fdHashMap[s2WordTable[framePos]]--) + 1;
}

std::string
D2::reverseComplement(std::string sequence) {
    std::reverse(sequence.begin(), sequence.end());
    const char Complements[] = {'T', 'C', 'G', 'A'};
    for (int i = 0; i < ((int) sequence.size()); i++) {
        sequence[i] = Complements[(int) CharToInt[(int) sequence[i]]];
    }
    return sequence;
}

float
D2::analyze(const int otherEST) {
    if (otherEST == refESTidx) {
        return 0; // distance to self will be 0
    } else if ((otherEST >= 0) && (otherEST < EST::getESTCount())) {
        // Check with the heuristic chain
        if ((chain == NULL) || (chain->shouldAnalyze(otherEST))) {
            int sed = 0;
            int minSed = frameSize*4; // won't be exceeded
            int s1FramePos = 0;
            // Get sequences
            EST *estS1 = EST::getEST(refESTidx);
            EST *estS2 = EST::getEST(otherEST);
            std::string sq1 = estS1->getSequence();
            ASSERT ( sq1.size() > 0 );
            std::string sq2 = estS2->getSequence();
            ASSERT ( sq2.size() > 0);

            // Initialize the word table for s2
            buildWordTable(s2WordTable, sq2);
 
            // Initialize the frequency differential hashmap for s1-s2
            buildFdHashMaps(&sed);

            minSed = sed;

            const int numFramesS1 = sq1.size() - frameSize;
            const int numFramesS2 = sq2.size() - frameSize;

            // Main d2 algorithm (from Zimmermann paper)
            while (s1FramePos < numFramesS1) {
                if (s1FramePos != 0) {
                    for (int i = 0; (i < frameShift && s1FramePos++ < numFramesS1); i++) {
                        refShiftUpdateFd(&sed, s1FramePos);
                        if (sed < minSed) {
                            alignmentMetric = s1FramePos;
                        minSed = sed;
                        if (minSed == 0) {
                            break;
                        }
                        }
                    }
                }
                for (int s2FramePos = 1; s2FramePos <= numFramesS2; s2FramePos++) {
                    rightShiftUpdateFd(&sed, s2FramePos);
                    if (sed < minSed)  {
                        minSed = sed;
                        alignmentMetric = s1FramePos - s2FramePos;
                        if (minSed==0) {
                            break;
                        }
                    }
                }
                if (s1FramePos != numFramesS1) {
                    for (int i = 0; (i < frameShift && s1FramePos++ < numFramesS1); i++) {
                        refShiftUpdateFd(&sed, s1FramePos);
                        if (sed < minSed) {
                            minSed = sed;
                            alignmentMetric = s1FramePos;
                            if (minSed==0) {
                                break;
                            }
                        }
                    }
                }
                for (int s2FramePos = numFramesS2-1; s2FramePos >= 0; s2FramePos--) {
                    leftShiftUpdateFd(&sed, s2FramePos);
                    if (sed < minSed) {
                        minSed = sed;
                        alignmentMetric = s1FramePos - s2FramePos;
                        if (minSed==0) {
                            break;
                        }
                    }
                }
            }

            //printf("%d %d %d\n", refESTidx, otherEST, minSed);
    
            return minSed;
        } else {
            // Heuristic chain said no, so return a dummy distance
            return (4*frameSize);
        }
    } else {
        // Invalid est index.
        return -1;
    }
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
