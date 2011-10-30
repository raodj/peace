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
#include "MPIStats.h"

#include <sstream>

// Some defines used as MPI Send/Recv tags in
// managerFormGlobalContig() and workerFormGlobalContig() methods
#define CONTIG_COUNTER_VALUES_TAG  200
#define NT_DISTRIBUTIONS_TAG       201
#define OVERHANGS_TAG              202

ContigMaker::ContigMaker(const int initialESTidx, const ESTList* cDnaList,
                         const bool tallyRoot)
    : rootESTidx(initialESTidx),
      rootSeq(cDnaList->get(initialESTidx)->getSequence()),
      rootSeqParentPos(0), estList(cDnaList) {
    ASSERT ( estList != NULL );
    ASSERT ( (rootESTidx >= 0) && (rootESTidx < estList->size()) );

    // Initialize logical bounds of the ntDistributions vector
    leftMostNtPos  = 0;
    rightMostNtPos = 0;
    rootESTNtPos   = 0;
    // Reset the overhangs
    leftOverhang        = "";
    leftOverhangESTidx  = -1;
    rightOverhang       = "";
    rightOverhangESTidx = -1;
    // Now add the root EST to our initial contig right at the center
    // of the contig with room to grow.
    rootESTLen = rootSeq.size();
    BatonAlignmentInfo rootInfo;
    rootInfo.setInfo(rootESTidx, -1, 0, false);
    rootInfo.getSegmentList().push_back(Segment(0, 0, rootESTLen));
    // Add the root info to our alignment information
    add(rootInfo, tallyRoot);
    // Add a suitable base segment entry as well.
    if (tallyRoot) {
        baseSegments.push_back(BaseSegmentInfo(rootESTidx, 0, rootESTLen - 1));
    }
}

ContigMaker::ContigMaker(const int baseSegIdx,
                         const std::string& rootESTSeq,
                         const int globalESTStartPos,
                         const int overhangStartPos,
                         const int overhangLen,
                         const ESTList* cDnaList,
                         const bool tallyRoot)
    : rootESTidx(-1),
      rootSeq(rootESTSeq),
      rootSeqParentPos(globalESTStartPos),
      estList(cDnaList) {
    ASSERT ( estList != NULL );
    // Initialize logical bounds of the ntDistributions vector
    rootESTLen     = rootSeq.size();    
    leftMostNtPos  = 0;
    rightMostNtPos = rootESTLen - 1;
    rootESTNtPos   = 0;
    // Reset the overhang nucleotide sequences
    leftOverhang        = "";
    leftOverhangESTidx  = -1;
    rightOverhang       = "";
    rightOverhangESTidx = -1;
    // Ensure sufficient capacity in ntDistributions for further growth.
    ensureCapacity(leftMostNtPos, rightMostNtPos);
    // Tally just the overhanging sub-fragment (as needed)
    if (tallyRoot) {
        tally(Segment(overhangStartPos,overhangStartPos, overhangLen),rootSeq);
        // Add a suitable base segment entry for it.
        baseSegments.push_back(BaseSegmentInfo(baseSegIdx, overhangStartPos,
                                       overhangStartPos + overhangLen - 1));
    }
}

ContigMaker::~ContigMaker() {
    // Nothing to be done for now.
}

void
ContigMaker::add(const BatonAlignmentInfo& ai) {
    // Sanity check to ensure we have at least one segment to work with
    ASSERT( ai.getSegmentList().size() > 0 );
    // Let overloaded method do the actual work
    add(ai, true);
}

void
ContigMaker::add(const BatonAlignmentInfo& ai, const bool tallyBases) {
    const EST* othEST = (ai.getESTIndex() != -1 ?
                         estList->get(ai.getESTIndex()) : NULL);
    // Update the left-most nucleotide position based on overhang
    const SegmentList &segList = ai.getSegmentList();
    const Segment &firstSeg    = segList[0];
    // Update the left-most nucleotide position
    leftMostNtPos = std::min(leftMostNtPos,
                             rootESTNtPos + firstSeg.getRefESTOffset());
    // Check and track left overhangs only when tallying bases
    if (tallyBases) {
        if ((firstSeg.getOthESTOffset() - firstSeg.getRefESTOffset()) >
            (int) leftOverhang.size()) {
            // We have a left-overhang. The new overhang has some
            // nucleotides to the left of the left-most overhang.  Track
            // only the longest left overhang.
            const std::string& othSeq =
                (othEST != NULL ? othEST->getSequence() : "");
            leftOverhang = othSeq.substr(0, firstSeg.getOthESTOffset() -
                                         firstSeg.getRefESTOffset());
            leftOverhangESTidx = (othEST != NULL ? othEST->getID() : -1);
        }
    }

    // Track the right edges including right-most nucleotides and overhang
    const Segment &lastSeg = segList[segList.size() - 1];
    // Track right-most nt
    rightMostNtPos = std::max(rightMostNtPos,
                              rootESTNtPos + lastSeg.getRelativeOffset() +
                              lastSeg.getLength() - 1);
    // Track the right overhang    
    if (tallyBases) {
        const int rightOverhangLen = othEST == NULL ? 0 :
            (othEST->getSequenceLength() - lastSeg.getOthESTEndOffset() - 1);
        if (rightOverhangLen > (int) rightOverhang.size()) {
            // We have a right-overhang. Track the nculeotides that don't
            // overlap.
            const std::string& othSeq =  othEST->getSequence();
            rightOverhang = othSeq.substr(lastSeg.getOthESTEndOffset() + 1);
            rightOverhangESTidx = (othEST != NULL ? othEST->getID() : -1);
        }
    }
    
    // Ensure we have sufficient capacity in ntDistributions to merge
    // the information for this alignment. In the call below the
    // leftMostNtPos, rightMostNtPos, rootESTNtPos, leftOverhangEnd,
    // rightOverhangStart get updated.
    ensureCapacity(leftMostNtPos, rightMostNtPos);
    
    // Track frequency of nucleotide occurrences using segment information
    // for non-reference fragments here.
    if (tallyBases) {
        // Get the regular or reverse-complement sequence depending on the
        // type of alignment.
        const std::string othSeq = (ai.isRCAlignment() ?
                                    othEST->getRCSequence() :
                                    othEST->getSequence());
        // Tally bases for each segment
        for(size_t i = 0; (i < segList.size()); i++) {
            tally(segList[i], othSeq);
        }
        // Save alignment information for future reference only when
        // real ESTs are involved even when we are tallying
        if (ai.getESTIndex() != -1) {
            alignedESTs.push_back(ai);
        }
    }
}

void
ContigMaker::tally(const Segment& seg, const std::string& othSeq) {
    const int othESTEndIdx = seg.getOthESTEndOffset();
    // Iterate over the nucleotides and track frequencies
    for(int othESTidx = seg.getOthESTOffset(),
            ntDestPos = rootESTNtPos + seg.getRefESTOffset();
        (othESTidx <= othESTEndIdx);
        othESTidx++, ntDestPos++) {
        ntDistributions[ntDestPos].add(othSeq[othESTidx]);
    }
}

void
ContigMaker::ensureCapacity(const int leftNtPos, const int rightNtPos,
                            const int padding) {
    if ((leftNtPos >= 0) && (rightNtPos < (int) ntDistributions.size())) {
        // We have sufficient capacity already
        return;
    }
    // We need to grow the ntDistributions vector to left and/or right
    // as needed.
    const NtDistr empty;
    // Pre-reserve memory in the vector for growth below
    const size_t FinalSize = ntDistributions.size() +
        (leftNtPos <= 0 ? (padding - leftNtPos) : 0) +
        (rightNtPos >= (int) ntDistributions.size() ?
         (padding + rightNtPos) : 0);
    ntDistributions.reserve(FinalSize);

    // First handle left-end as necessary
    if (leftNtPos <= 0) {
        // Shift entries to the right with padding entries open for
        // further growth to left
        const int rightShift = padding - leftNtPos;
        ntDistributions.insert(ntDistributions.begin(), rightShift, empty);
        // Adjust the starting positions of the various variables
        leftMostNtPos      += rightShift;
        rightMostNtPos     += rightShift;
        rootESTNtPos       += rightShift;
    }

    // Next handle right-end as necessary
    while (ntDistributions.size() < FinalSize) {
        ntDistributions.push_back(empty);
    }
}

int
ContigMaker::formGlobalContig() {
    const int MyRank       = MPI_GET_RANK();
    const int ProcessCount = MPI_GET_SIZE();
    
    if (MyRank != MANAGER_RANK) {
        // Perform worker tasks -- First send local contig information
        // to the manager.
        sendContigData(MANAGER_RANK, alignedESTs.size());
        // Receive global contig from the manager (don't merge but
        // replace)
        return recvContigData(MANAGER_RANK, false);
    }
    // When control drops here we are performing manager-process
    // operations. First gather local-contig information from all the
    // workers into this contig to form the global contig.
    int alignedESTsCount = alignedESTs.size();
    for(int procID = 1; (procID < ProcessCount); procID++) {
        // Receive data and merge it into this contig.
        alignedESTsCount += recvContigData(procID, true);
    }
    // Now we have global contig at manager. Distribute it to all the
    // workers.
    for(int procID = 1; (procID < ProcessCount); procID++) {
        // Send global contig data to a worker.
        sendContigData(procID, alignedESTsCount);
    }
    // Return total numebr of aligned ESTs.
    return alignedESTsCount;
}

void
ContigMaker::sendContigData(const int procID, const int alignedESTsCount) {
    // Form combined left+right overhang to send overhang nucleotides
    // in one shot.
    const std::string overhangs = leftOverhang + rightOverhang;
    // First send out various counters. The counters include
    // additional information to ease receiving ntDistributions and
    // overhang nucleotide sequences.
    const int counterValues[9] = {leftMostNtPos, rightMostNtPos,
                                  rootESTNtPos, leftOverhang.size(),
                                  overhangs.size(), alignedESTsCount,
                                  ntDistributions.size(), leftOverhangESTidx,
                                  rightOverhangESTidx};
    MPI_SEND(counterValues, 9, MPI_INT, procID, CONTIG_COUNTER_VALUES_TAG);
    // Next send the overhang nucleotide sequences
    MPI_SEND(&overhangs[0], overhangs.size(), MPI_CHAR, procID, OVERHANGS_TAG);
    // Finally send the ntDistributions array
    const int bytesToSend = ntDistributions.size() * sizeof(NtDistr);
    MPI_SEND(&ntDistributions[0], bytesToSend, MPI_CHAR,
             procID, NT_DISTRIBUTIONS_TAG);
}

int
ContigMaker::recvContigData(const int procID, const bool mergeData) {
    // The sending sequence in sendContigData() and the sequence of
    // receive operations in this method must be correctly
    // coordinated.

    // First receive the counter values.
    int counterValues[9];
    TRACK_IDLE_TIME(MPI_RECV(counterValues, 9, MPI_INT, procID,
                             CONTIG_COUNTER_VALUES_TAG));
    // Next receive the left+right overhang nucleotide
    const int& ovrHangSize = counterValues[4];
    std::string overhangs(ovrHangSize, '-');
    TRACK_IDLE_TIME(MPI_RECV(&overhangs[0], ovrHangSize, MPI_CHAR,
                             procID, OVERHANGS_TAG));
    // Next receive the nucleotide distribution data. If we are
    // merging we receive data in a separate array to minimize memory
    // operations.
    const int& ntCount = counterValues[6];
    if (mergeData) {
        NtDistrList otherNtDistributions(ntCount);
        TRACK_IDLE_TIME(MPI_RECV(&otherNtDistributions[0], ntCount * sizeof(NtDistr),
                                 MPI_CHAR, procID, NT_DISTRIBUTIONS_TAG));
        // Now merge the received data with our local information
        merge(otherNtDistributions, counterValues[0],
              counterValues[1], counterValues[2],
              overhangs.substr(0, counterValues[3]),
              counterValues[7],  // otherLeftOverhangESTidx
              overhangs.substr(counterValues[3]),
              counterValues[8]); // otherRightOverhangESTidx
    } else {
        // Overwrite situation (getting complete data from manager)
        ntDistributions.resize(ntCount);
        TRACK_IDLE_TIME(MPI_RECV(&ntDistributions[0], ntCount*sizeof(NtDistr),
                                 MPI_CHAR, procID, NT_DISTRIBUTIONS_TAG));
        // Update all counter values as well.
        leftMostNtPos       = counterValues[0];
        rightMostNtPos      = counterValues[1];
        rootESTNtPos        = counterValues[2];
        leftOverhang        = overhangs.substr(0, counterValues[3]);
        rightOverhang       = overhangs.substr(counterValues[3]);
        leftOverhangESTidx  = counterValues[7];
        rightOverhangESTidx = counterValues[8];
    }
    // Return the number of aligned ESTs reported by the sender.
    return counterValues[5];
}

std::string
ContigMaker::getConsensus() const {
    return getConsensus(leftMostNtPos, rightMostNtPos);
}

std::string
ContigMaker::getNtDistribution(const std::string& separator) const {
    std::ostringstream result;
    for(int pos = leftMostNtPos; (pos <= rightMostNtPos); pos++) {
        const int* freq     = ntDistributions[pos].freq;
        const int  freqBase = std::max_element(freq, freq + 4) - freq;
        result << (pos > leftMostNtPos ? separator : "")
               << freq[freqBase];
    }
    return result.str();
}

std::string
ContigMaker::getConsensus(const int startPos, const int endPos) const {
    ASSERT((startPos >= leftMostNtPos) && (startPos <= rightMostNtPos));
    ASSERT((endPos   >= leftMostNtPos) && (endPos   <= rightMostNtPos));
    ASSERT(endPos    >= startPos);
    // Form the consensus sequence for the requested range
    std::string consensus;
    consensus.reserve(endPos - startPos + 1);
    for(int i = startPos; (i <= endPos); i++) {
        const char nt = ntDistributions[i].getConsensus();
        consensus.push_back(nt);
    }
    // Consensues sequence for the range is ready!
    return consensus;
}

void
ContigMaker::merge(const NtDistrList& otherNtDistributions,
                   const int otherLeftMostNtPos, const int otherRightMostNtPos,
                   const int otherRootESTNtPos,
                   const std::string& otherLeftOverhang,
                   const int otherLeftOverhangESTidx,
                   const std::string& otherRightOverhang,
                   const int otherRightOverhangESTidx) {
    // Compute difference in root EST position (which serves as our
    // anchor position for merging).
    const int deltaRootESTNtPos = rootESTNtPos - otherRootESTNtPos;
    // Update the left and right overhang information prior to merge
    // (while we still have our current, unmodified counters). In this
    // method read "adj" as "adjsuted"
    const int adjOtherLeftMostNtPos = otherLeftMostNtPos + deltaRootESTNtPos;
    // Check and track the left-most nucleotide sequences
    if (adjOtherLeftMostNtPos < leftMostNtPos) {
        // The other contig has nucleotides all the way to the
        // left. So simply retain its left-overhang as the final
        // overhang.
        leftOverhang = otherLeftOverhang;
        leftOverhangESTidx = otherLeftOverhangESTidx;
    } else if (adjOtherLeftMostNtPos == leftMostNtPos) {
        // Both have the same left-most position. Just choose the
        // longest overhang.
        if (leftOverhang.size() < otherLeftOverhang.size()) {
            leftOverhang       = otherLeftOverhang;
            leftOverhangESTidx = otherLeftOverhangESTidx;
        }
    }
    // Appropriately handle the right overhang information
    const int adjOtherRightMostNtPos = otherRightMostNtPos + deltaRootESTNtPos;
    if (adjOtherRightMostNtPos > rightMostNtPos) {
        // The other contig is more right than this one. So simply
        // retain its right-overhang as the final overhang.
        rightOverhang       = otherRightOverhang;
        rightOverhangESTidx = otherRightOverhangESTidx;
    } else if (adjOtherRightMostNtPos == rightMostNtPos) {
        // Both have the same right-most position. Just choose the
        // longest overhang.
        if (rightOverhang.size() < otherRightOverhang.size()) {
            rightOverhang       = otherRightOverhang;
            rightOverhangESTidx = otherRightOverhangESTidx;
        }
    }
    
    // Ensure we have sufficient capacity to merge the
    // information. The following method call updates: leftMostNtPos,
    // rightMostNtPos, rootESTNtPos
    ensureCapacity(adjOtherLeftMostNtPos, adjOtherRightMostNtPos);
    // Next check and update our left-most and right-most nucleotide
    // positions. We have to recompute delta root position as
    // rootESTNtPos may have changed in call to ensureCapacity()
    const int adjDeltaPos = rootESTNtPos - otherRootESTNtPos;
    leftMostNtPos  = std::min(leftMostNtPos,
                              otherLeftMostNtPos + adjDeltaPos);
    rightMostNtPos = std::max(rightMostNtPos,
                              otherRightMostNtPos + adjDeltaPos);
    // Since the above operations may have updated rootESTNtPos, we
    // need to compute adjusted values.
    const int adjOtherStartIndex = rootESTNtPos - otherRootESTNtPos;
    for(int i = otherLeftMostNtPos; (i <= otherRightMostNtPos); i++) {
        ntDistributions[i + adjOtherStartIndex].add(otherNtDistributions[i]);
    }
}

ContigMaker
ContigMaker::getOverhangContigMaker(const bool useLeftOverhang,
                                    const int  windowSize,
                                    const bool tallyRoot) const {
    // The following two variables will eventually contain the root
    // nucleotide sequence and its starting position with which we are
    // going to construct the new contig.
    std::string ovrHngSeq;  // The overhang nucleotide sequence
    int ovrHngSeqPos = -1;  // The position of overhang sequence
    // Common variable used below
    const int overhangLen    = (useLeftOverhang ? leftOverhang.size() :
                                rightOverhang.size());
    const int overhangESTidx = (useLeftOverhang ? leftOverhangESTidx :
                                rightOverhangESTidx);
    // The following logic is applied to both left and right overhang
    // below:
    //
    //    1. If overhang len is shorter than 1 window size add
    //       sufficient consensus nucleotides to make ovrHngSeq
    //       to be 1 window long.
    //
    //    2. Otherwise add 1/2 window long consensus nucleotides to
    //       the overhang.
    if (useLeftOverhang) {
        // Left overhang case (consensus sequence to be obtained from
        // left-end of this contig and added to the right of the
        // overhang).
        // Compute consensus end position.
        const int suffixEndPos = leftMostNtPos +
            (overhangLen >= windowSize ?  (windowSize / 2) :
             (windowSize - overhangLen)) - 1;
        // Get consensus suffix (can be an empty string).
        const std::string suffix =
            getConsensus(leftMostNtPos, std::min(rightMostNtPos, suffixEndPos));
        // Update the variables common to left-or-right case.
        ovrHngSeqPos = (leftMostNtPos - overhangLen);
        ovrHngSeq    = leftOverhang  + suffix;
    } else {
        // Right overhang case (consensus sequence to be obtained from
        // right-end of this contig and added to the left of the
        // overhang).
        const int overhangLen  = rightOverhang.size();
        // Compute consensus start position
        const int prefixStartPos = rightMostNtPos -
            (overhangLen >= windowSize ? (windowSize / 2) - 1:
             (windowSize - overhangLen));
        // Get the consensus sequence prefix (can be empty string)
        const std::string prefix =
            getConsensus(std::max(leftMostNtPos, prefixStartPos),
                         rightMostNtPos);
        // Update the variables common to left-or-right case.
        ovrHngSeq    = prefix + rightOverhang;
        ovrHngSeqPos = rightMostNtPos - prefix.size() + 1;
    }
    ASSERT ( ovrHngSeq.size() > 0 );
    // Create a new ContigMaker with the necessary information.
    const int overhangStartPos = (useLeftOverhang ? 0 :
                                  (ovrHngSeq.size() - overhangLen));
    ContigMaker cm(overhangESTidx, ovrHngSeq, ovrHngSeqPos, overhangStartPos,
                   overhangLen, estList, tallyRoot);
    return cm;
}

void
ContigMaker::merge(const ContigMaker& subContig) {
    // Merge alignment information with segment offsets suitably
    // adjusted.
    alignedESTs.reserve(alignedESTs.size() + subContig.alignedESTs.size());
    for(size_t i = 0; (i < subContig.alignedESTs.size()); i++) {
        if (subContig.alignedESTs[i].getESTIndex() != -1) {
            alignedESTs.push_back(subContig.alignedESTs[i]);
            alignedESTs.back().adjustRefOffsets(subContig.rootSeqParentPos -
                                                rootESTNtPos);
        }
    }
    // Next merge nucleotides. The merge happens such that the
    // nulceotides at subContig.rootESTNtPos align with this (parent)
    // contig at position subContig.rootSeqParentPos
    const int deltaPos = rootESTNtPos - subContig.rootSeqParentPos;
    // merge the nucleotides from the subContig into this (parent) contig.
    merge(subContig.ntDistributions, subContig.leftMostNtPos,
          subContig.rightMostNtPos,  subContig.rootESTNtPos + deltaPos,
          subContig.leftOverhang,    subContig.leftOverhangESTidx,
          subContig.rightOverhang,   subContig.rightOverhangESTidx);
    // Merge the base segment information from the sub-contig to
    // this main contig.
    baseSegments.reserve(baseSegments.size() + subContig.baseSegments.size());
    for(size_t i = 0; (i < subContig.baseSegments.size()); i++) {
        const BaseSegmentInfo& bsi = subContig.baseSegments[i];
        baseSegments.push_back(BaseSegmentInfo(bsi.getESTIndex(),
                                               bsi.getStartPos() - deltaPos,
                                               bsi.getEndPos()   - deltaPos));
    }
}

void
ContigMaker::prettyPrint(std::ostream& os) const {
    // First print the reference fragment at the appropriate location
    if (rootESTidx == -1) {
        os << "Ref Seq   :" << std::string(rootESTNtPos - leftMostNtPos, ' ')
           << rootSeq      << std::endl;
    }
    // Print each alignment at the appropriate location
    for(size_t i = 0; (i < alignedESTs.size()); i++) {
        if (alignedESTs[i].getESTIndex() == -1) {
            // Dummy EST for root sequence. Skip it.
            continue;
        }
        const int estIdx = alignedESTs[i].getESTIndex();
        char estInfo[16];
        sprintf(estInfo, "EST#%-6d:", estIdx);
        os << estInfo;
        const std::string seq = estList->get(estIdx)->getSequence();
        alignedESTs[i].prettyPrint(os, seq, rootESTNtPos - leftMostNtPos);
        os << std::endl;
    }
    // Finally print the consensus sequence for the contig.
    std::string consensus = getConsensus(leftMostNtPos, rightMostNtPos);
    std::cout << "Consensus :" << consensus << std::endl;
}

void
ContigMaker::populateContig(Contig& contig, const bool setConsensus,
                            const bool setNtDistributions) const {
    if (setConsensus) {
        contig.setConsensus(getConsensus());
    }
    if (setNtDistributions) {
        contig.setNtDistribution(ntDistributions.begin(),
                                 ntDistributions.begin());
    }
    // Next add alignment information to the contig. We need to adjust
    // the left position such that the nucleotide at the left-most
    // position is at logical position zero.
    const int adjustedLeftPos = rootESTNtPos - leftMostNtPos;
    for(size_t i = 0; (i < alignedESTs.size()); i++) {
        const BatonAlignmentInfo& bai = alignedESTs[i];
        const EST* est                = estList->get(bai.getESTIndex());
        // Collect the required alignment information.
        const std::string cigar = bai.getCIGAR(est->getSequence());
        // Update start position so that the left-most nucleotide will
        // occupy logical position of zero.
        const int position = bai.getStartOffset() + adjustedLeftPos;
        ASSERT ( position >= 0 );
        // Create the generic alignment object from baton-alignment
        AlignmentInfo ai(bai.getESTIndex(), bai.isRCAlignment(),
                         position, cigar, 255, est->getSequenceLength(),
                         bai.getScore());  // score is optional
        // Add it to the contig.
        contig.addAlignmentInfo(ai);
    }
    // Finally, add the base segment information to the
    // contig. Similar to the alignment, we need to adjust the left
    // position such that the nucleotide at the left-most position is
    // at logical position zero (for this contig).
    for(size_t i = 0; (i < baseSegments.size()); i++) {
        // Get current base segment information.
        const BaseSegmentInfo& bsi = baseSegments[i];
        // Create new adjusted base segment information
        BaseSegmentInfo adjustedBSI(bsi.getESTIndex(),
                                    bsi.getStartPos() + adjustedLeftPos,
                                    bsi.getEndPos()   + adjustedLeftPos);
        // Add the base segment information to the contig
        contig.addBaseSegmentInfo(adjustedBSI);
    }
}

#endif
