#ifndef BATON_LIST_CPP
#define BATON_LIST_CPP

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

#include "BatonList.h"
#include "ESTCodec.h"

#include <algorithm>
#include <iterator>
#include <cstdio>

// The static bitShift value. This is currently fixed to 2 as each
// nucleotide is encoded using a 2-bit encoding.
int BatonList::bitShift = 4;

// The bit-mask to be used to drop bits that are no longer needed.
// This value is set based on the n-mers used for baton heads.
int BatonList::bitMask = 0;

BatonList::BatonList(const EST* srcEst, const int nMers, const bool makeRC,
                     const int windowSize) :
    est(srcEst),  nMerSize(nMers), isRC(makeRC), sequence("") {
    // Now that the constant instance variables have been initialized,
    // let's get the helper method to actually build the batons for us.
    buildBatons(windowSize);
}

BatonList::BatonList(const char *seq, const int nMers, const bool makeRC,
                     const int windowSize)
    : est(NULL), nMerSize(nMers), isRC(makeRC), sequence(seq) {
    // Now that the constant instance variables have been initialized,
    // let's get the helper method to actually build the batons for us.
    buildBatons(windowSize);    
}

void
BatonList::getCodedNMers(IntVector& codedNMers) const {
    // First, setup the bit mask based on the number of
    // n-mers. Compute bit mask that will retain only the bits
    // corresponding to a given word size.  Each entry in a word takes
    // up 2 bits and that is why the following formula involves a 2.
    bitMask = (1 << (nMerSize * 2)) - 1;
    // Now, build n-mer list using normal or reverse complement encoder
    if (isRC) {
        ESTCodec::RevCompEncoder<bitShift, bitMask> encoder;
        getCodedNMers(getSequence(), codedNMers, encoder);
    } else {
        ESTCodec::NormalEncoder<bitShift, bitMask> encoder;
        getCodedNMers(getSequence(), codedNMers, encoder);
    }    
}

void
BatonList::buildBatons(const int winSize) {
    // Ensure we have at least one window.
    const int seqLen = getSequenceLength() - nMerSize + 1;
    windowCount      = seqLen / winSize; // number of windows (can be 0!)
    const int excess = seqLen % winSize;  // bases that don't evenly fit
    // If more than 50% of windowSize does not fit then create an
    // additional window by decreasing window size.  Otherwise, grow
    // window size by a few bases to evenly distribute the excess.
    if (excess > winSize / 2) {
        // Create a new window and make window size smaller.
        windowCount++;  // can't be 0 anymore.
        windowSize = seqLen / windowCount;
    } else if (windowCount > 0) {
        // Only a few bases don't seem to fit. Grow each window by a
        // few bases to accommodate the excess.
        windowSize = winSize + (excess / windowCount);
    } else if (windowCount == 0) {
        ASSERT( seqLen < winSize );
        windowSize = seqLen;
    }
    
    // Due to overlapping windows we actually have twice the number of
    // windows, assumming the sequence is long enough.
    windowCount = std::max(1, windowCount * 2);
    // Next get the encoded n-mers for the given sequence via helper
    // method. The helper method handles both normal and reverse
    // complement cases.
    IntVector codedNMers;
    getCodedNMers(codedNMers);
    // Set up the lastOccurrence array (to track last occurence of a
    // given n-mer) with all possible n-mer combinations.  The last
    // occurrence values initialized to -1.
    const int MaxNmerValue = (1 << (nMerSize * 2));
    IntVector lastOccurrence(MaxNmerValue, -1);
    // Initialize the batons list to hold batons whose baton heads
    // correspond to various n-mer values.
    batons.resize(MaxNmerValue);
    // Now compute the nMerCode for each n-mer and create batons
    const int MerCount = codedNMers.size();
    for(int idx = 0; (idx < MerCount); idx++) {
        const int nMerCode     = codedNMers[idx];
        const int lastOccurIdx = lastOccurrence[nMerCode];
        if (lastOccurIdx != -1) {
            // Found a new baton whose head corresponds to the
            // nMerCode value.
            batons[nMerCode].push_back(Baton(lastOccurIdx,
                                             idx - lastOccurIdx + nMerSize));
        }
        // Track the index where we noticed the curren n-mer to
        // aid baton during future iterations of this for-loop
        lastOccurrence[nMerCode] = idx;
    }
    // Now sort the list of batons in each entry in the batons list
    // based on the baton lengths to optimize baton-based comparisons
    for(size_t i = 0; (i < batons.size()); i++) {
        // Sort based on length.  The Baton::operator<() is internally
        // used by std::sort for comparisons.
        std::stable_sort(batons[i].begin(), batons[i].end());
    }
}

void
BatonList::getWindowPairs(const BatonList& other,
                          std::vector<WindowPair>& pairList,
                          const int threshold) const {
    // Create the temporary vector to hold baton frequency counts
    std::vector< IntVector > identicalBatonCount;
    // Let the other method do the actual computations and populate
    // identicalBatonCount 2-D array.
    getWindowPairs(other, identicalBatonCount, threshold);
    // Now process the entries in the 2-D array and create WindowPair
    // objects of entries that exceed the given threshold.
    for(size_t win1 = 0; (win1 < identicalBatonCount.size()); win1++) {
        const IntVector& win1Entries = identicalBatonCount[win1];
        for(size_t win2 = 0; (win2 < win1Entries.size()); win2++) {
            if (win1Entries[win2] >= threshold) {
                // Found a pair of windows that have a sufficient
                // number of identical batons between them.
                pairList.push_back(WindowPair(win1, win2, win1Entries[win2]));
            }
        }
    }
}

void
BatonList::getWindowPairs(const BatonList& other,
                          std::priority_queue<WindowPair>& pairList,
                          const int threshold) const {
    // Create the temporary vector to hold baton frequency counts
    std::vector< IntVector > identicalBatonCount;
    // Let the other method do the actual computations and populate
    // identicalBatonCount 2-D array.
    getWindowPairs(other, identicalBatonCount, threshold);
    // Now process the entries in the 2-D array and create WindowPair
    // objects of entries that exceed the given threshold.
    for(size_t win1 = 0; (win1 < identicalBatonCount.size()); win1++) {
        const IntVector& win1Entries = identicalBatonCount[win1];
        for(size_t win2 = 0; (win2 < win1Entries.size()); win2++) {
            if (win1Entries[win2] >= threshold) {
                // Found a pair of windows that have a sufficient
                // number of identical batons between them.
                pairList.push(WindowPair(win1, win2, win1Entries[win2]));
            }
        }
    }
}

// Method to identify window-pairs whose identical baton count exceeds
// the given threshold value.
void
BatonList::getWindowPairs(const BatonList& other,
                          std::vector< IntVector >& identicalBatonCount,
                          const int threshold) const {
    const int winSize1 = windowSize / 2;  // This is half width
    const int winSize2 = other.windowSize / 2;
    ASSERT ( threshold > 0 );
    
    // The following 2-d array (implemented using vector-of-vectors)
    // is used to track the number of identical batons between pairs
    // of overlapping windows that we can seen so far.  Note the +1 in
    // windowCount -- This is needed because we divide winSize1 and
    // winSize2 by two.  This causes the last few nucleotides to
    // overflow window count. So we conservatively add 1 to the number
    // of elements in the matrix identicalBatonCount.
    const IntVector ZeroVector(other.windowCount + 1, 0);
    identicalBatonCount.clear();   // Clear out old entries
    identicalBatonCount.resize(windowCount + 1, ZeroVector); // Create 2-d vector
    // For each n-mer code, compare baton lengths to find matching batons
    for(size_t nMerCode = 0; (nMerCode < batons.size()); nMerCode++) {
        // Note that here the batons are sorted based on their
        // lengths. So the search for matching batons with the same
        // n-mer ends proceeds differently.
        NmerBatonList::const_iterator baton1 = this->batons[nMerCode].begin();
        NmerBatonList::const_iterator baton2 = other.batons[nMerCode].begin();
        // Search for matching batons until we exhaust a list
        while ((baton1 != this->batons[nMerCode].end()) &&
               (baton2 != other.batons[nMerCode].end())) {
            if (baton1->getLength() == baton2->getLength()) {
                // Found matching batons! Update identical baton
                // counts and add section pair to pairList if count
                // exceeds the specified threshold value.
                tallyBatons(baton1->getStartIndex() / winSize1, windowCount,
                            baton2->getStartIndex() / winSize2,
                            other.windowCount, identicalBatonCount);
                // Onto the next pair of batons to check
                baton1++;
                baton2++;
            } else if (baton1->getLength() < baton2->getLength()) {
                // Baton2 is longer. So a matching baton should be
                // further down in this->batons[nMerCode] list.
                baton1++;
            } else {
                // Baton1 is longer. So a matching baton should be
                // further down in other.batons[nMerCode] list.                
                baton2++;
            }
        }
    }
}

void
BatonList::tallyBatons(const int baton1Win, const int numWin1,
                       const int baton2Win, const int numWin2,
                       std::vector< IntVector > &identicalBatonCount) const {
    // The following checks are broken down into the following three
    // main categories: (i) one of the batons are in the first (or
    // left-most) window; (ii) one of the batons are in the last (or
    // right-most) window; (iii) both batons are somewhere in the
    // middle (that is, neither the first nor the last).
    if ((baton1Win == 0) || (baton2Win == 0)) {
        // Case (i): One or both batons are in the first window
        if (baton1Win == baton2Win) { // Both are zero!
            identicalBatonCount[0][0]++;
        } else if (baton1Win == 0) {
            identicalBatonCount[0][baton2Win]++;
            identicalBatonCount[0][baton2Win - 1]++;            
        } else {
            ASSERT( baton2Win == 0 );
            identicalBatonCount[baton1Win][0]++;
            identicalBatonCount[baton1Win - 1][0]++;
        }
    } else if ((baton1Win >= numWin1) || (baton2Win >= numWin2)) {
        // Case (ii): One or both batons are in the last window
        if ((baton1Win >= numWin1) && (baton2Win >= numWin2)) {
            identicalBatonCount[baton1Win - 1][baton2Win - 1]++;
        } else if (baton1Win >= numWin1) {
            identicalBatonCount[baton1Win - 1][baton2Win]++;
            identicalBatonCount[baton1Win - 1][baton2Win - 1]++;
        } else {
            identicalBatonCount[baton1Win][baton2Win - 1]++;
            identicalBatonCount[baton1Win - 1][baton2Win - 1]++;
        }
    } else {
        // Case (iii): both batons are somewhere in the middle
        identicalBatonCount[baton1Win - 1][baton2Win - 1]++;
    }
}

std::ostream&
operator<<(std::ostream& os, const BatonList& batonList) {
    // Print the batons in each n-mer list.
    for(size_t i = 0; (i < batonList.batons.size()); i++) {
        // First print the baton head n-mer sequence.
        char buffer[4];
        sprintf(buffer, "%02x: ", (unsigned int) i);
        os << buffer;
        // Print each baton in the baton list.
        const NmerBatonList& nMerList = batonList.batons[i];
        std::copy(nMerList.begin(), nMerList.end(),
                  std::ostream_iterator<Baton>(os, " "));
        // End the line of the baton.
        os << std::endl;
    }
    // Return reference to the stream as per API requirement
    return os;
}

#endif
