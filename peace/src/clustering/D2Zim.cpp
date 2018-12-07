#ifndef D2_ZIM_CPP
#define D2_ZIM_CPP

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

#include "D2Zim.h"
#include "ESTList.h"
#include "ESTCodec.h"
#include "ArgParser.h"
#include "HeuristicChain.h"
#include <algorithm>

// Since the following variables are passed in a template parameter,
// it is defined to be static (to ensure that it has external linkage
// as per the ISO/ANSI standards requirement).
int D2Zim::BitMask = 0;
// A copy of the wordSize instance variable (in the base class) used
// for template parameter
int D2Zim::statWordSize = 0;

D2Zim::D2Zim() : FWAnalyzer("D2") {
    fdHashMap     = NULL;
    s1WordTable   = NULL;
    s2WordTable   = NULL;

    frameShift    = 1;
    wordTableSize = 0;
    NumWordsWin   = 0;    
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
D2Zim::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add its parameters first    
    FWAnalyzer::addCommandLineArguments(argParser);
    // The set of arguments for this class.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--frameShift", "Frame Shift (default=1)",
         &frameShift, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Setup valid arguments for this analyzer
    argParser.addValidArguments(ArgsList);
}

bool
D2Zim::initialize() {
    if (!FWAnalyzer::initialize()) {
        // Error occured when initalizing base class.  This is no good.
        return false;
    }
    // Compute the size of the frequency differential hashmap
    const int MapSize = 1 << (wordSize * 2);

    // Compute bit mask that will retain only the bits corresponding
    // to a given word size.  Each entry in a word takes up 2 bits and
    // that is why the following formula involves a 2.
    BitMask = (1 << (wordSize * 2)) - 1;
    // Setup NumWordsWin instance variable
    NumWordsWin = frameSize - wordSize;
    
    // Initialize the fd hash map
    fdHashMap = new int[MapSize];

    // Compute word table size and initialize word tables
    ASSERT( estList != NULL );
    wordTableSize = estList->getMaxESTLen() + frameSize;
    s1WordTable   = new int[wordTableSize];
    s2WordTable   = new int[wordTableSize];
    
    // Everything went on well.
    return true;
}

int
D2Zim::setReferenceEST(const EST* est) {
    // Call corresponding method in heuristic chain
    if (chain != NULL) {
        chain->setReferenceEST(est);
    }
    if (est != NULL) {
        refEST = est;
        // init ref-est word table
        const char* s1 = est->getSequence();
        buildWordTable(s1WordTable, s1);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

void
D2Zim::buildWordTable(int* wordTable, const char* s) {
    // The following word size is a static instance so that it can be
    // passed as template parameter

    // init ref-est word table
    memset(wordTable, 0, sizeof(int) * wordTableSize);
    // First compute the hash for a single word.
    ASSERT ( s != NULL );
    const int length = strlen(s) - wordSize;
    statWordSize = wordSize;
    ESTCodec::NormalEncoder<statWordSize, BitMask> encoder;
    register int hash = 0;
    int ignoreMask    = 0;
    for(int i = 0; (i < wordSize); i++) {
        hash = encoder(hash, s[i], ignoreMask);
    }
    // Setup the word table entry for the first word.
    wordTable[0] = hash;
    // Fill in the word table
    for (int i = 1; (i <= length); i++) {
        hash = encoder(hash, s[i], ignoreMask);
        if (!ignoreMask) {
            // If ignoreMask is zero that indicates that the hash is
            // not tainted due to 'n' bp in sequence and it is good.
            wordTable[i] = hash;
        }
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
D2Zim::getMetric(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    if (otherEST->getID() == refEST->getID()) {
        return 0; // distance to self will be 0
    }
    // Ensure other EST is fully populated
    estList->repopulate(otherEST);
    
    // Perform operations for D2
    int sed = 0;
    int minSed = frameSize * 4; // won't be exceeded
    int s1FramePos = 0;
    // Get sequences
    const char* sq1  = refEST->getSequence();
    ASSERT ( sq1 != NULL );
    const char* sq2 = otherEST->getSequence();
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
    
    return (float) minSed;
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
