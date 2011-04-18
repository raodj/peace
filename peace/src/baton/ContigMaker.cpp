#ifndef CONTIG_MAKER_CPP
#define CONTIG_MAKER_CPP

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

#include "ContigMaker.h"
#include "ESTCodec.h"

#include <algorithm>

ContigMaker::ContigMaker(const int rootESTidx,
                         const ESTList* list) : estList(list) {
    ASSERT( estList != NULL );
    // Save the root for future reference
    AlignmentInfo root(rootESTidx, rootESTidx, 0, 0, false);
    alignedESTs.push_back(root);
    // Setup other instance variables
    currRefEntry       = -1;
    leftMostEntry      = 0;  // Index in alignedESTs vector
    leftMostProcessed  = false;
    rightMostEntry     = 0; // Index in alignedESTs vector
    rightMostProcessed = true;  // Already processed
    leftOrRight        = false; // Start with left.
}

ContigMaker::~ContigMaker() {
    // Nothing to be done in the destructor for now.
}

int
ContigMaker::nextReference(int& entry, bool& processed) {
    if ((!processed) && (entry != -1) && (entry < (int) alignedESTs.size())) {
        processed  = true;  // Set flag to indcate entry is processed
        return entry;       // Return the entry back
    }
    // If this entry has already been processed and/or it is invalid
    // we return -1 to indicate an invalid reference EST index.
    return -1;
}

int
ContigMaker::nextReference() {
    // First choice giving full precedence to right-most (if
    // leftOrRight is true) or left-most entry.
    currRefEntry = (leftOrRight ? nextReference(rightMostEntry,
                                                rightMostProcessed) :
                    nextReference(leftMostEntry, leftMostProcessed));
    // If the above default did not find an entry, then try the
    // opposite options.
    if (currRefEntry == -1) {
        currRefEntry = (!leftOrRight ? nextReference(rightMostEntry,
                                                     rightMostProcessed) :
                        nextReference(leftMostEntry, leftMostProcessed));
    }
    // Switch the default direction to check the next time this method
    // is called.
    leftOrRight = !leftOrRight;
    // Return the current reference entry we have. This value may be
    // -1 if both left and right have been processed and no new
    // entries exceeding this range have been added.
    return (currRefEntry != -1) ? alignedESTs[currRefEntry].getESTIndex() : -1;
}

void
ContigMaker::add(const AlignmentInfo& info) {
    const int AlignedESTCount = (int) alignedESTs.size();
    ASSERT ( (currRefEntry >= 0) && (currRefEntry < AlignedESTCount) );
    ASSERT ( currRefEntry != -1 );
    // Create new entry for the est being added.
    alignedESTs.push_back(info);
    // Update the left-most or right-most entry based on the offset.
    ASSERT ((leftMostEntry >= 0) && (leftMostEntry < AlignedESTCount));
    if (alignedESTs[leftMostEntry].getAlignmentPos() > info.getAlignmentPos()){
        // The newly added fragment is more left than the previous one
        leftMostEntry     = alignedESTs.size() - 1;
        leftMostProcessed = false;
    }
    // Compute the logical offset occupied by the right-most base
    // to decide if the newly added entry has the right-most one.
    ASSERT ((rightMostEntry >= 0) && (rightMostEntry < AlignedESTCount));
    if (alignedESTs[rightMostEntry].getRightMostNtPos(estList) <
        info.getRightMostNtPos(estList)) {
        // The newly added entry has nucleotides to the right of
        // the earlier right-most entry. So set the new one as the
        // right-most entry.
        rightMostEntry     = alignedESTs.size() - 1;
        // Set right most processed to false only if the left and
        // right entry are different. If they are the same we process
        // only the left most entry once.
        rightMostProcessed = (leftMostEntry == rightMostEntry);
    }
}

void
ContigMaker::formContig(const bool verbose) {
    if (alignedESTs.empty()) {
        // Nothing to form a gene out of.  This happens when this
        // method gets called again from the destructor to ensure all
        // the contigs are actually flushed out.
        return;
    }
    // Obtain the complete range of offsets for the contig to be formed.
    const int LeftMostOffset  = alignedESTs[leftMostEntry].getAlignmentPos();
    const int RightMostOffset =
        alignedESTs[rightMostEntry].getRightMostNtPos(estList);
    // Print some information if in verbose mode.
    if (verbose) {
        std::cout << "Forming contig using: "  << alignedESTs.size()
                  << " genes in the range: "   << LeftMostOffset
                  << " to "                    << RightMostOffset
                  << std::endl;
    }
    // Create a two dimensional array to count the number of bases
    // occuring at each logical position.
    TwoDArray baseDistribution;
    // Get base distributions populated via a helper method to
    // streamline the code.
    getBaseDistributions(baseDistribution, verbose);
    // Now that we have the base distributions for each position,
    // process the distributions to choose the nucleotide with the
    // highest occurrence count to form the consensus sequence.
    std::string consensus = formConsensus(baseDistribution);
    // Print the consensus sequence
    std::cout << "Consensus sequence for contig:\n"
              << consensus << std::endl;
}

void
ContigMaker::getBaseDistributions(TwoDArray& distributions,
                                     const bool verbose) const {
    // Obtain the complete range of offsets for the contig to be
    // formed. Get the left-most (LM) and right-most (RM) offsets
    const int LMOffset = alignedESTs[leftMostEntry].getAlignmentPos();
    const int RMOffset = alignedESTs[rightMostEntry].getRightMostNtPos(estList);

    // First populate the distributions array with all zeros.
    const std::vector<int> ZeroVector(5, 0);  // Five for ATCG and N
    // Create 2-d vector
    distributions.resize(RMOffset - LMOffset + 1, ZeroVector);
    // Compute the base distributions using the aligned fragments
    for(size_t index = 0; (index < alignedESTs.size()); index++) {
        // Get the sequence and sequence length for convenience.
        const EST  *est = estList->get(alignedESTs[index].getESTIndex());
        std::string seq = est->getSequence();
        if (alignedESTs[index].isRCAlignment()) {
            seq = makeRevComp(seq);
        }

        // Determine the actual starting position for this fragment
        // within the nucleotide position of the contig being formed.
        const int startPos = alignedESTs[index].getAlignmentPos() - LMOffset;
        if (verbose) {
            std::cout << "EST #" << index << "\t" << std::string(startPos, ' ')
                      << seq << std::endl;
        }
        // Update the base distribution..
        const int seqLen = seq.size();
        for(int seqIdx = 0; (seqIdx < seqLen); seqIdx++) {
            distributions[startPos + seqIdx][ESTCodec::encode(seq[seqIdx])]++;
        }
    }
}

std::string
ContigMaker::makeRevComp(const std::string& seq) const {
    // Create a reversed-copy of the incoming string.
    std::string rc(seq.rbegin(), seq.rend());
    // Flip each nucleotide in the sequence to its complement.
    for(size_t idx = 0; (idx < rc.size()); idx++) {
        rc[idx] = ESTCodec::getComp(rc[idx]);
    }
    // Return the reverse complement version.
    return rc;
}

std::string
ContigMaker::formConsensus(const TwoDArray& distributions) const {
    // Create consensus sequence with a series of '-'s
    std::string consensus(distributions.size(), '-');
    // For each nucleotide index (ntIdx), determine the nucleotide
    // that has the highest frequency of occurrence. Use that value to
    // determine the base pair character at the given position.
    for(size_t ntIdx = 0; (ntIdx < distributions.size()); ntIdx++) {
        // Determine the index of the element with the highest occurrence.
        const int codedMaxBase = std::max_element(distributions[ntIdx].begin(),
                                                  distributions[ntIdx].end()) -
            distributions[ntIdx].begin();
        // Now codedMaxBase has value in the range 0 to 4
        // corresponding to ATCG-.  Convert the value into a base
        // character
        consensus[ntIdx] = ESTCodec::decode(codedMaxBase);
    }
    // Return the consensus sequence to the caller
    return consensus;
}

#endif
