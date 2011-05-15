#ifndef DEFAULT_SEQUENCE_ALIGNER_CPP
#define DEFAULT_SEQUENCE_ALIGNER_CPP

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

#include "DefaultSequenceAligner.h"
#include "BatonListCache.h"
#include "AlignmentInfo.h"
#include <limits>

DefaultSequenceAligner::DefaultSequenceAligner() :
    SequenceAligner("default") {
    // Initialize command-line arguments to default values
    threshold          = 7;  // aka (criticalValue)
    numPermittedErrs   = 10;
    goodAlignmentScore = 15;
    // Clear out the reference sequence baton information
    refEST       = NULL;
    refBatonList = NULL;
}

DefaultSequenceAligner::~DefaultSequenceAligner() {
}

void
DefaultSequenceAligner::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    SequenceAligner::addCommandLineArguments(argParser);
    // Now add our custom parameters.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--threshold", "Min. baton count above which 2 reads are similar",
         &threshold, ArgParser::INTEGER},        
        {"--permErrs", "#of differences to be accepted when checking candidate alignments",
         &numPermittedErrs, ArgParser::INTEGER},        
        {"--goodScore", "A good-enough alignment score to short circuit exhaustive alignment analysis",
         &goodAlignmentScore, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(ArgsList);
}

bool
DefaultSequenceAligner::initialize() {
    // Let the base class do its initialization
    if (!SequenceAligner::initialize()) {
        // Error occured when initializing.  This is no good.
        return false;
    }
    // Ensure threshold is valid.
    if (threshold < 1) {
        std::cerr << getName()
                  << ": Threshold must be >= 1 (suggested value: 7)"
                  << "(use --threshold option)\n";
        return false;
    }
    // Check and update the goodalignmentScore
    if (goodAlignmentScore == -1) {
        goodAlignmentScore = std::numeric_limits<int>::max();
    }    
    // Initialization successful.
    return true;
}

int
DefaultSequenceAligner::setReferenceEST(const EST* est) {
    // Clear out any baton list created for a direct-given consensus
    // type sequence.
    if (refBatonList != NULL) {
        // delete refBatonList;
        refBatonList = NULL;
        refEST       = NULL;
    }
    // Setup reference EST for future use
    refEST = est;
    if (est != NULL) {
        // Pre-compute normal baton list, if it is not already in our
        // blCache...
        refBatonList = blCache->getBatonList(est, false);
        ASSERT ( refBatonList != NULL );
        // Build up the n-mer list for future use.
        refBatonList->getCodedNMers(codedRefNmers);
        return 0; // everything went well
    }
    // Invalid EST index.
    return 1;
}

// Primary method method to compute the alignment between the reference
// sequence and a given EST.
bool
DefaultSequenceAligner::align(const EST* otherEST, AlignmentInfo& info) {
    // We should have a reference baton list already set
    ASSERT( refBatonList != NULL );
    // The following variables are populated in the align() method
    // calls below.
    int refAlignPos, othAlignPos, score;
    // First check to see if alignment with regular (not
    // reverse complement) nucleotide sequence of otherEST yeilds good
    // results.
    if (align(otherEST, refAlignPos, othAlignPos, score, false)) {
        // Found a sufficiently good alignment. Update info and return.
        info = AlignmentInfo(otherEST->getID(), refEST->getID(),
                             refAlignPos - othAlignPos, score, false);
        return true;
    }
    // Since first round of analysis failed, next check to see if
    // alignment with reverse complement) nucleotide sequence of
    // otherEST yeilds good results.
    if (align(otherEST, refAlignPos, othAlignPos, score, true)) {
        // Found a sufficiently good alignment. Update info and return
        info = AlignmentInfo(otherEST->getID(), refEST->getID(),
                             refAlignPos - othAlignPos, score, true);
        return true;
    }
    // Both normal and reverse complement analysis failed.
    return false;
}

// Helper method method to compute the alignment between the reference
// sequence and a given EST by searching window-pairs with a large
// number of identical batons
bool
DefaultSequenceAligner::align(const EST* otherEST, int& refAlignPos,
                              int& othAlignPos, int& score, const bool doRC) {
    // We should have a reference baton list already set
    ASSERT( refBatonList != NULL );
    // Build the normal or reverse complement (based on doRC flag)
    // baton list for the otherEST for analysis.
    const BatonList *othBatonList = blCache->getBatonList(otherEST, doRC);
    ASSERT ( othBatonList != NULL );
    // check to see if we have sufficient number of identical batons
    // within a logical window in the two sequences.
    std::vector<WindowPair> pairList;
    refBatonList->getWindowPairs(*othBatonList, pairList, threshold);
    if (pairList.empty()) {
        // No good window pairs were found. No need to try to align
        // any further in this case
        return false;
    }
    // Obtain nucleotide sequences for reference below
    const std::string refSeq(refEST->getSequence());
    const std::string otherSeq(otherEST->getSequence());
    
    // Set starting default score to zero
    score = 0;
    // Search among pairs of windows (with sufficent number of
    // identical batons) returned in the list for best alignment.  But
    // first get the coded n-mers for the otherEST.
    IntVector codedOthNmers;
    othBatonList->getCodedNMers(codedOthNmers);
    for(size_t winIdx = 0; (winIdx < pairList.size()); winIdx++) {
        // Get alignment based on the current window-pair
        int currRefAlignPos, currOthAlignPos;
        int currScore =
            getBestAlignmentInWindow(refSeq, *refBatonList, codedRefNmers,
                                     pairList[winIdx].getWindow1(),
                                     currRefAlignPos, otherSeq,
                                     *othBatonList, codedOthNmers,
                                     pairList[winIdx].getWindow2(),
                                     currOthAlignPos,
                                     numPermittedErrs, goodAlignmentScore);
        // Track the best alignment we have had so far...
        if (currScore > score) {
            score       = currScore;
            refAlignPos = currRefAlignPos;
            othAlignPos = currOthAlignPos;
        }
        // Short circuit exhaustive alignment search if a reasonably
        // good alignment was found. The parameter goodAlignmentScore is
        // also-known-as signficiantAlignment.
        if (score >= goodAlignmentScore) {
            break;
        }
    }
    // Indicate we have an alignment
    return true;
}

// Helper method that determines best alignment using batons in a
// given window-pair (called only by align() method above).
int
DefaultSequenceAligner::getBestAlignmentInWindow(const std::string& refSeq,
                                        const BatonList& refBatonList,
                                        const IntVector& codedRefNmers,
                                        const int refWindow,
                                        int&  refAlignPos,
                                        const std::string& otherSeq,
                                        const BatonList& othBatonList,
                                        const IntVector& codedOthNmers,
                                        const int othWindow,
                                        int&  othAlignPos,
                                        const int numPermittedErrs,
                                        const int goodAlignmentScore) const {
    // Verify consistency
    ASSERT ( refBatonList.getNmerSize() == othBatonList.getNmerSize() );
    // Compute all possible encoding values for the n-mers
    const int MaxNmerValue = (1 << (refBatonList.getNmerSize() * 2));
    // For each baton pair that falls within the respective windows in
    // the two baton lists, track the best alignment thus far..
    int bestScore = 0;
    for(int nMer = 0; (nMer < MaxNmerValue); nMer++) {
        for(NmerBatonList::const_iterator refBaton = refBatonList[nMer].begin();
            (refBaton != refBatonList[nMer].end()); refBaton++) {
            if (!refBatonList.isBatonInWindow(*refBaton, refWindow)) {
                // This baton is not in the window we are
                // checking. So skip it.
                continue;
            }
            for(NmerBatonList::const_iterator othBaton = othBatonList[nMer].begin();
                (othBaton != othBatonList[nMer].end()); othBaton++) {
                if (!othBatonList.isBatonInWindow(*othBaton, othWindow)) {
                    // This baton is not in the window we are
                    // checking. So skip it.
                    continue;
                }
                // Check to ensure that the batons are identical.
                if (refBaton->getLength() != othBaton->getLength()) {
                    // These two batons have the same k-mer heads but
                    // are of different lengths. Consequently, they
                    // are not identical batons. Skip this pair for
                    // further analysis.
                    continue;
                }
                // Determine alignment score for the current pair of
                // refBaton and othBaton
                int currScore =
                    getAlignmentScore(refSeq, codedRefNmers,
                                      refBaton->getStartIndex(), otherSeq,
                                      codedOthNmers, othBaton->getStartIndex(),
                                      numPermittedErrs);
                // Track best score and update baton start positions.
                if (currScore > bestScore) {
                    // Found  a better score. So update it.
                    bestScore   = currScore;
                    refAlignPos = refBaton->getStartIndex();
                    othAlignPos = othBaton->getStartIndex();
                }
                // If we fould a good enough alignment then exit right
                // away without exploring other (possibly better)
                // alignments, favoring run time..
                if (bestScore >= goodAlignmentScore) {
                    return bestScore;
                }
            }
        }
    } // n-mer for-loop
    // Return the best alignment we found
    return bestScore;
}

// This is a custom helper method that is used to search to the
// left of an identical baton-pair starting positions to locate
// first position where two consecutive nucleotide differences
// occur.
Segment
DefaultSequenceAligner::makeSegment(const std::string& refSeq,
                                    const int baton1Pos,
                                    const int refSeqStartPos,
                                    const int refSeqEndPos,
                                    const std::string& otherSeq,
                                    const int baton2Pos,
                                    const int othSeqStartPos,
                                    const int othSeqEndPos) const {
    // Use the helper methods to search to the left to determine the
    // left-bounds of Segment
    const int leftSpanLen =
        countMatchingBasesToLeft(refSeq,   baton1Pos, refSeqStartPos,
                                 otherSeq, baton2Pos, othSeqStartPos);
    // Use the helper methods to search to the right to determine the
    // right-bounds of Segment
    const int rightSpanLen =
        countMatchingBasesToRight(refSeq,   baton1Pos, refSeqEndPos,
                                  otherSeq, baton2Pos, othSeqEndPos);
    // Now form the segment using the left and right spans.
    return Segment(baton1Pos - leftSpanLen, baton2Pos - rightSpanLen,
                   leftSpanLen + rightSpanLen);
}

// This method attempts to figure out the left-position where the
// nucleotides in refSeq and otherSeq differ by two consecutively
// different base-pairs.
int
DefaultSequenceAligner::countMatchingBasesToLeft(const std::string& refSeq,
                                                 const int baton1Pos,
                                                 const int refSeqStartPos,
                                                 const std::string& otherSeq,
                                                 const int baton2Pos,
                                                 const int othSeqStartPos) const {
    // Ensure we can even move to the left if possible.
    if ((baton1Pos <= refSeqStartPos) || (baton2Pos <= othSeqStartPos)) {
        // The baton positions are already at or before start of
        // sequences
        return 0;
    }
    // Setup variables used in the while loop below to determine
    // number of different nucleotides seen thus far.
    register int refSeqPos = baton1Pos - 1;  // index in refSeq
    register int othSeqPos = baton1Pos - 1;  // index int otherSeq
    register int diffNtCnt = 0;              // Count of different nts thus far
    
    // This loop typically ends when either currRefPos or currOthPos
    // becomes -1 (zero is a valid index).
    while ((refSeqPos >= refSeqStartPos) && (othSeqPos >= othSeqStartPos)) {
        if (refSeq[refSeqPos] == otherSeq[othSeqPos]) {
            // Nucleotides at this position match. Reset diffNtCnt
            // variable to indicate there are no differences thusfar.
            diffNtCnt = 0;
        } else {
            // Found different nucleotides at this position.
            diffNtCnt++;
            if (diffNtCnt > 1) {
                // When the diffNtCnt is greater than 1, then we found
                // two consecutive nucleotide differences. Time to
                // break out of the while loop.
                break;
            }
        }
        // Check the previous nucleotides for (mis)match
        refSeqPos--;
        othSeqPos--;
    }
    // Return number of positions moved to left in reference sequence
    // (in case you are wondering, it will be the same in other
    // sequence as well). However, we never return a negative value.
    return baton1Pos - refSeqPos;
}

// This method attempts to figure out the right-position where the
// nucleotides in refSeq and otherSeq differ by two consecutively
// different base-pairs.
int
DefaultSequenceAligner::countMatchingBasesToRight(const std::string& refSeq,
                                                  const int baton1Pos,
                                                  const int refSeqEndPos,
                                                  const std::string& otherSeq,
                                                  const int baton2Pos,
                                                  const int othSeqEndPos) const {
    // Ensure we can even move to the right if possible.
    if ((baton1Pos >= refSeqEndPos) || (baton2Pos >= othSeqEndPos)) {
        // The baton positions are already at or before start of
        // sequences
        return 0;
    }
    // Setup variables used in the while loop below to determine
    // number of different nucleotides seen thus far.
    register int refSeqPos = baton1Pos + 1;  // index in refSeq
    register int othSeqPos = baton1Pos + 1;  // index int otherSeq
    register int diffNtCnt = 0;              // Count of different nts thus far
    
    // This loop typically ends when either currRefPos or currOthPos
    // becomes greater than refSeqEndPos or othSeqEndPos respectively.
    while ((refSeqPos < refSeqEndPos) && (othSeqPos < othSeqEndPos)) {
        if (refSeq[refSeqPos] == otherSeq[othSeqPos]) {
            // Nucleotides at this position match. Reset diffNtCnt
            // variable to indicate there are no differences thusfar.
            diffNtCnt = 0;
        } else {
            // Found different nucleotides at this position.
            diffNtCnt++;
            if (diffNtCnt > 1) {
                // When the diffNtCnt is greater than 1, then we found
                // two consecutive nucleotide differences. Time to
                // break out of the while loop.
                break;
            }
        }
        // Check the next pair of nucleotides for (mis)match
        refSeqPos++;
        othSeqPos++;
    }
    // Return number of positions moved to left in reference sequence
    // (in case you are wondering, it will be the same in other
    // sequence as well). However, we never return a negative value.
    return refSeqPos - baton1Pos;
}

// Helper method to determine start (or seed) for a segment to the
// left of a given pair (ref for reference sequence, and oth for other
// sequence being aligned with ref) of positions in two cDNA
// fragments.
bool
DefaultSequenceAligner::findSeedToLeft(const IntVector& codedRefNmers,
                                       const int refStartPos,
                                       const IntVector& codedOthNmers,
                                       const int othStartPos,
                                       int& refSeedPos, int& othSeedPos) const {
    // Compute the left-position where we are going to stop searching
    const int refLeftStopPos = (refStartPos > 10) ? (refStartPos - 9) :
        refStartPos;
    const int othLeftStopPos = (othStartPos > 10) ? (othStartPos - 9) :
        refStartPos;
    // Check to see if we have something to search for
    if ((refLeftStopPos <= 0) || (othLeftStopPos <= 0)) {
        // No seeds to search for.
        return false;
    }
    // 
}


int
DefaultSequenceAligner::getAlignmentScore(const std::string& UNUSED(refSeq),
                                          const IntVector& codedRefNmers,
                                          int baton1Pos,
                                          const std::string& UNUSED(otherSeq),
                                          const IntVector& codedOthNmers,
                                          int baton2Pos,
                                          int numPermittedErrs) const {
    // Set the maximum length of alignment we can achieve based on the
    // shorter of the two cDNA fragments.
    const int MaxAlignLen = std::min(codedRefNmers.size(),
                                     codedOthNmers.size());
    int leftMoves  = 0; // count of bases we have managed to move left  so far
    int rightMoves = 0; // count of bases we have managed to move right so far
    
    // The following two flags track if we would like to further
    // extend/expore alignment towards the left or towards the right
    bool tryMoveToLeft  = true;
    bool tryMoveToRight = true;

    // The following loop runs until we exhaust the quota of permitted
    // number of errors (aka differences) or we achieve the maximum
    // alignment possible.
    do {
        // Check how far we can go either left and/or right until we
        // encounter a difference.
        const int rightStepsTaken = (!tryMoveToRight) ? 0 :
            tryRightExtension(codedRefNmers, baton1Pos + rightMoves,
                              codedOthNmers, baton2Pos + rightMoves);
        const int leftStepsTaken  = (!tryMoveToLeft) ? 0 :
            tryLeftExtension(codedRefNmers, baton1Pos - leftMoves,
                             codedOthNmers, baton2Pos + leftMoves);
        // Check and explore the direction in which we moved a lot
        // more further.
        if (rightStepsTaken >= leftStepsTaken) {
            // Track how far right, across iterations we have moved so far
            rightMoves    += rightStepsTaken;
            tryMoveToLeft  = false;     // Don't go left in next iteration
            tryMoveToRight = true;      // Try going right more next
        } else {
            // Track how far left, across iteration we have moved so far
            leftMoves     += leftStepsTaken; 
            tryMoveToLeft  = true;     // Try going left more in next iteration
            tryMoveToRight = false;    // Don't go right in next iteration
        }
        // When control drops here we have encountered one difference
        // either to the left or to the right. So track it.
        numPermittedErrs--;
    } while ((numPermittedErrs > 0) && ((leftMoves + rightMoves)< MaxAlignLen));
    // Return the sum of the left and right distances we have managed
    // to move while accruing numPermittedErrs number of differences/errors.
    return (leftMoves + rightMoves);
}

int
DefaultSequenceAligner::tryLeftExtension(const IntVector& codedRefNmers,
                                         int refIndexPos,
                                         const IntVector& codedOthNmers,
                                         int othIndexPos) const {
    // Save starting point to compute delta-moves for return value below.
    const int startRefIndexPos = refIndexPos; 
    while ((refIndexPos > 0) && (othIndexPos > 0)) {
        refIndexPos--; // Move one step left in reference sequence
        othIndexPos--; // Move one step left in other
        if (codedRefNmers[refIndexPos] != codedOthNmers[othIndexPos]) {
            break;
        }
    }
    return startRefIndexPos - refIndexPos;
}

int
DefaultSequenceAligner::tryRightExtension(const IntVector& codedRefNmers,
                                          size_t refIndexPos,
                                          const IntVector& codedOthNmers,
                                          size_t othIndexPos) const {
    if ((refIndexPos >= codedRefNmers.size()) ||
        (othIndexPos >= codedOthNmers.size())) {
        // One of the indexes is too large. Cannot scan to right.
        return 0;
    }
    // Save starting point to compute delta-moves for return value below.
    const size_t startRefIndexPos = refIndexPos;
    do {
        refIndexPos++; // Move one step right in reference sequence
        othIndexPos++; // Move one step right in other
    } while ((refIndexPos < codedRefNmers.size()) &&
             (othIndexPos < codedOthNmers.size()) &&
             (codedRefNmers[refIndexPos] != codedOthNmers[othIndexPos]));
    // Return the delta moved to the right. It will at least be one.
    return refIndexPos - startRefIndexPos;
}

#endif

//  LocalWords:  signficiant
