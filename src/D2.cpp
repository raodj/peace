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
//
//---------------------------------------------------------------------------

#include "D2.h"
#include "EST.h"
#include <cmath>

// The simple translation table to convert chars to int.
char D2::CharToInt[256];

int D2::MapSize = (int) pow(4, wordSize);

int* fdHashMap;
int* rcFdHashMap;


// Compute bit mask that will retain only the bits corresponding
// to a given word size.  Each entry in a word takes up 2 bits and
// that is why the following formula involves a 2.
int D2::BitMask = (1 << (wordSize * 2)) - 1;

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
    // Initialize the fd hash maps
    // These are hash maps that map hashes of words to integers, just like
    // in CLU.  However, in this case, we are keeping track of the frequency
    // difference.  So fdHashMap contains the frequency differences for each
    // word in s1 and s2; rcFdHashMap contains the frequency differences for
    // each word in s1 and the RC of s2.  This is done according to the D2
    // algorithm pseudocode in Zimmermann's paper.
    fdHashMap = new int[MapSize];
    rcFdHashMap = new int[MapSize];
}

D2::~D2() {
    if (fdHashMap != NULL) {
        delete [] fdHashMap;
    }
    if (rcFdHashMap != NULL) {
        delete [] rcFdHashMap;
    }
}

void
D2::showArguments(std::ostream& os) {
    FWAnalyzer::showArguments(os);
    // Currently this class does not have any specific command line
    // arguments.  Consequently, there is nothing further to be done.
}

bool
D2::parseArguments(int& argc, char **argv) {
    // Currently this class does not have any specific command line
    // arguments.  The base class does all the necessary tasks.

    // Let the base class do processing and return the result.
    return FWAnalyzer::parseArguments(argc, argv);
}

int
D2::initialize() {
    int result = FWAnalyzer::initialize();
    if (result != 0) {
        // Error occured when loading ESTs.  This is no good.
        return result;
    }
    // Everything went on well.
    return 0;
}

int
D2::setReferenceEST(const int estIdx) {
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        return 0;
    }
    // Invalid est index.
    return 1;
}

void
D2::buildFdHashMaps(const char* sequence, int* sed, int* rcSed, int* leftHash,
		    int* rightHash) {
    // Clear out any old entries in the hash maps.
    memset(fdHashMap, 0, sizeof(int) * MapSize);
    memset(rcFdHashMap, 0, sizeof(int) * MapSize);
    // First compute the hash for a single word.
    ASSERT ( sequence != NULL );
    int hash = 0;
    int i;
    for(i = 0; (i < wordSize); i++) {
        hash <<= 2;
        hash  += CharToInt[(int) sequence[i]];
    }
    // Setup the hash entry for the first word.
    fdHashMap[hash]++;
    rcFdHashMap[hash]++;
    (*sed)++;
    (*rcSed)++;
    *leftHash = hash;

    for (; (i < frameSize); i++) { // only go to frame size
        hash <<= 2;
        hash  += CharToInt[(int) sequence[i]];
        hash  &= BitMask;
	// update sed and fd for this word
	*sed+=(fdHashMap[hash]++)*2+1;
	*rcSed+=(rcFdHashMap[hash]++)*2+1;
    }
    *rightHash = hash;
}

void
D2::initialUpdateFd(const char* sequence, int* sed, int* leftHash,
		    int* rightHash, bool rc) {
  // First compute the hash for a single word.
  ASSERT ( sequence != NULL );
  int hash = 0;
  int i;
  for(i = 0; (i < wordSize); i++) {
    hash <<= 2;
    hash  += CharToInt[(int) sequence[i]];
  }
  // update sed and fd for the first word
  rc ? *sed+=((rcFdHashMap[hash]--)*-2)+1 : *sed+=((fdHashMap[hash]--)*-2)+1;
  *leftHash = hash;

  for (; (i < frameSize); i++) { // only go to frame size
    hash <<= 2;
    hash  += CharToInt[(int) sequence[i]];
    hash  &= BitMask;
    // update sed and fd for this word
    rc ? *sed+=((rcFdHashMap[hash]--)*-2)+1 : *sed+=((fdHashMap[hash]--)*-2)+1;
  }
  *rightHash = hash;
}

void
D2::refShiftUpdateFd(const char* sequence, int* sed, int* rcSed,
		     int* leftHash, int* rightHash, const int framePos) {
  // update sed and fd from leftmost word falling out
  *sed+=((fdHashMap[*leftHash]--)*-2)+1;
  *rcSed+=((rcFdHashMap[*leftHash]--)*-2)+1;

  // update leftmost word's hash
  *leftHash <<= 2;
  *leftHash += CharToInt[(int) sequence[framePos+(wordSize-1)]];
  *leftHash &= BitMask;

  // update rightmost word's hash (new rightmost word)
  *rightHash <<= 2;
  //ASSERT(framePos+frameSize < strlen(sequence));
  *rightHash += CharToInt[(int) sequence[framePos+(frameSize-1)]];
  *rightHash &= BitMask;

  // update sed and fd from new rightmost word
  *sed+=((fdHashMap[*rightHash]++)*2)+1;
  *rcSed+=((rcFdHashMap[*rightHash]++)*2)+1;
}

void
D2::rShiftUpdateFd(const char* sequence, int* sed, int* leftHash, 
		   int* rightHash, const int framePos, bool rc) {
  // update sed and fd from leftmost word falling out
  rc ? *sed+=((rcFdHashMap[*leftHash]++)*2)+1 : 
    *sed+=((fdHashMap[*leftHash]++)*2)+1;

  // update leftmost word's hash
  *leftHash <<= 2;
  *leftHash += CharToInt[(int) sequence[framePos+(wordSize-1)]];
  *leftHash &= BitMask;

  // update rightmost word's hash (new rightmost word)
  *rightHash <<= 2;
  //ASSERT(framePos+frameSize < strlen(sequence));
  *rightHash += CharToInt[(int) sequence[framePos+(frameSize-1)]];
  *rightHash &= BitMask;

  // update sed and fd from new rightmost word
  rc ? *sed+=((rcFdHashMap[*rightHash]--)*-2)+1 : 
    *sed+=((fdHashMap[*rightHash]--)*-2)+1;
}

void
D2::lShiftUpdateFd(const char* sequence, int* sed, int* leftHash, 
		   int* rightHash, const int framePos, bool rc) {
  // update sed and fd from rightmost word falling out
  rc ? *sed+=((rcFdHashMap[*rightHash]++)*2)+1 : 
    *sed+=((fdHashMap[*rightHash]++)*2)+1;

  // update rightmost word's hash
  *rightHash >>= 2;
  int temp = CharToInt[(int) sequence[framePos+frameSize-wordSize]];
  temp <<= 2 * (wordSize-1);
  *rightHash += temp;
  *rightHash &= BitMask;

  // update leftmost word's hash (new leftmost word)
  *leftHash >>= 2;
  temp = CharToInt[(int) sequence[framePos]];
  temp <<= 2 * (wordSize-1);
  *leftHash += temp;
  *leftHash &= BitMask;

  // update sed and fd from new leftmost word
  rc ? *sed+=((rcFdHashMap[*leftHash]--)*-2)+1 : 
    *sed+=((fdHashMap[*leftHash]--)*-2)+1;
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
    int sed = 0; int rcSed = 0;
    int minSed = 200; int rcMinSed = 200;
    int s1FramePos = 0;
    int s1FrameLeftHash = 0; int s1FrameRightHash = 0;
    int s2FrameLeftHash = 0; int s2FrameRightHash = 0;
    int rcFrameLeftHash = 0; int rcFrameRightHash = 0;
    // Get sequences
    EST *estS1 = EST::getEST(refESTidx);
    EST *estS2 = EST::getEST(otherEST);
    std::string sequence = estS1->getSequence();
    ASSERT ( sequence.size() > 0 );
    const char* sq1 = sequence.c_str();
    std::string sequence2=estS2->getSequence();
    ASSERT(sequence2.size() > 0);
    const char* sq2 = sequence2.c_str();

    // create RC of S2
    const char* sqRC = reverseComplement(sequence2).c_str();

    //printf("%d %d\n", wordSize, frameSize);
 
    // Initialize the frequency differential hashmaps for s1-s2 and s1-rcS2
    buildFdHashMaps(sq1, &sed, &rcSed, &s1FrameLeftHash, &s1FrameRightHash);
    //printf("%d %d\n", sed, rcSed);
    // Update using the second sequence in the pairs of sequences
    initialUpdateFd(sq2, &sed, &s2FrameLeftHash, &s2FrameRightHash, false);
    initialUpdateFd(sqRC, &rcSed, &rcFrameLeftHash, &rcFrameRightHash, true);

    //printf("%d %d\n", sed, rcSed);

    const int numFramesS1 = strlen(sq1) - frameSize;
    const int numFramesS2 = strlen(sq2) - frameSize;
    // numWordsRC_S2 is same as numWordsS2

    // Main d2 algorithm (from Zimmermann paper)
    while (s1FramePos < numFramesS1) {
      if (s1FramePos != 0) {
	s1FramePos++;
	refShiftUpdateFd(sq1, &sed, &rcSed, &s1FrameLeftHash,
			 &s1FrameRightHash, s1FramePos);
      }
      for (int s2FramePos = 1; s2FramePos <= numFramesS2; s2FramePos++) {
	rShiftUpdateFd(sq2, &sed, &s2FrameLeftHash, &s2FrameRightHash,
		      s2FramePos, false);
        rShiftUpdateFd(sqRC, &rcSed, &rcFrameLeftHash, &rcFrameRightHash, 
		      s2FramePos, true);
	if (sed < minSed) minSed = sed;
	if (rcSed < rcMinSed) rcMinSed = rcSed;
      }
      if (s1FramePos != numFramesS1) {
	s1FramePos++;
	refShiftUpdateFd(sq1, &sed, &rcSed, &s1FrameLeftHash, 
			 &s1FrameRightHash, s1FramePos);
      }
      for (int s2FramePos = numFramesS2-1; s2FramePos >= 0; s2FramePos--) {
	lShiftUpdateFd(sq2, &sed, &s2FrameLeftHash, &s2FrameRightHash,
		      s2FramePos, false);
	lShiftUpdateFd(sqRC, &rcSed, &rcFrameLeftHash, &rcFrameRightHash, 
		      s2FramePos, true);
	if (sed < minSed) minSed = sed;
	if (rcSed < rcMinSed) rcMinSed = rcSed;
      }
    }
    if (rcMinSed < minSed) {
      //printf("%d %d %d\n", refESTidx, otherEST, rcMinSed);
      return rcMinSed;
    } else {
      // printf("%d %d %d\n", refESTidx, otherEST, minSed);
      return minSed;
    }
       
  } else {
    // Invalid est index.
    return -1;
  }
}

//int
//D2::analyze() {
//    return -1;
//}
#endif
