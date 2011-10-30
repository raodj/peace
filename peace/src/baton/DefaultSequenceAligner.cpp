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
#include "BatonAlignmentInfo.h"

#include <algorithm>
#include <utility>
#include <limits>

DefaultSequenceAligner::DefaultSequenceAligner() :
    SequenceAligner("default") {
    // Initialize command-line arguments to default values
    threshold          = 7;  // aka (criticalValue)
    numPermittedErrs   = 10;
    goodAlignmentScore = 25; // aka significantAlignmentNumber
    minAlignScore      = 15; // Minimum alignment score needed before a
    // a alignment is declared as being good.
    // Clear out the reference sequence baton information
    refEST       = NULL;
    refBatonList = NULL;
}

DefaultSequenceAligner::~DefaultSequenceAligner() {
    if (refBatonList != NULL) {
        delete refBatonList;
        refBatonList = NULL;
        refEST       = NULL;
    }    
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
        {"--goodScore", "A good-enough seed-alignment score to short circuit exhaustive alignment analysis",
         &goodAlignmentScore, ArgParser::INTEGER},
        {"--minAlignScore", "A good-enough all-segments-alignment score to short circuit exhaustive alignment analysis",
         &minAlignScore, ArgParser::INTEGER},
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
    refSeq = refEST->getSequence();
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

int
DefaultSequenceAligner::setReferenceEST(const std::string& refcDNA) {
    // Clear out any baton list created for a direct-given consensus
    // type sequence.
    if (refBatonList != NULL) {
        if (refEST == NULL) {
            // This is a baton list created using a directly-specified
            // reference sequence.
            delete refBatonList;
        }
        refBatonList = NULL;
        refEST       = NULL;
    }
    // Setup reference sequence and baton list for future use
    refSeq = refcDNA;
    // Pre-compute normal baton list for reference sequence.
    refBatonList = new BatonList(refSeq.c_str(),
                                 blCache->getBatonHeadSize(),
                                 false,  // normal (not reference complement)
                                 blCache->getWindowSize());
    
    ASSERT ( refBatonList != NULL );
    // Build up the n-mer list for future use.
    refBatonList->getCodedNMers(codedRefNmers);
    return 0; // everything went well
}

// Primary method method to compute the alignment between the reference
// sequence and a given EST.
bool
DefaultSequenceAligner::align(const EST* otherEST,
                              BatonAlignmentInfo& alignInfo) {
    // We should have a reference baton list already set
    ASSERT( refBatonList != NULL );
    // A constant to make the code a bit more clear
    const bool doRevComp = true;
    // First check to see if alignment with regular (not
    // reverse complement) nucleotide sequence of otherEST yeilds good
    // results.
    if (align(otherEST, !doRevComp, alignInfo)) {
        // Found a sufficiently good alignment. Update info and return.
        return true;
    }
    // Since first round of analysis failed, next check to see if
    // alignment with reverse complement) nucleotide sequence of
    // otherEST yeilds good results.
    if (align(otherEST, doRevComp, alignInfo)) {
        return true;
    }
    // Both normal and reverse complement analysis failed.
    return false;
}

// Helper method method to compute the alignment between the reference
// sequence and a given EST by searching window-pairs with a large
// number of identical batons
bool
DefaultSequenceAligner::align(const EST* otherEST, const bool doRC,
                              BatonAlignmentInfo& alignInfo) {
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
    const std::string otherSeq(otherEST->getSequence());
    
    // Search among pairs of windows (with sufficent number of
    // identical batons) returned in the list for best alignment.  But
    // first get the coded n-mers for the otherEST.
    IntVector codedOthNmers;
    othBatonList->getCodedNMers(codedOthNmers);
    for(size_t winIdx = 0; (winIdx < pairList.size()); winIdx++) {
        // Get alignment based on the current window-pair
        int score =
            getBestAlignmentInWindow(refSeq, *refBatonList, codedRefNmers,
                                     pairList[winIdx].getWindow1(),
                                     otherSeq, *othBatonList, codedOthNmers,
                                     pairList[winIdx].getWindow2(),
                                     numPermittedErrs, goodAlignmentScore,
                                     alignInfo.getSegmentList());
        // Short circuit exhaustive alignment search if a reasonably
        // good alignment was found. The parameter goodAlignmentScore is
        // also-known-as significiantAlignment.
        if (score >= minAlignScore) {
            // Found an alignment. Update alignment information.
            alignInfo.setInfo(otherEST->getID(),
                              (refEST != NULL ? refEST->getID() : -1),
                              score, doRC);
            // Indicate we have an alignment
            return true;
        }
    }
    // Indicate we did not find a good-enough alignment
    return false;
}

// Helper method that determines best alignment using batons in a
// given window-pair (called only by align() method above).
int
DefaultSequenceAligner::getBestAlignmentInWindow(const std::string& refSeq,
                                        const BatonList& refBatonList,
                                        const IntVector& codedRefNmers,
                                        const int refWindow,
                                        const std::string& otherSeq,
                                        const BatonList& othBatonList,
                                        const IntVector& codedOthNmers,
                                        const int othWindow,
                                        const int UNUSED(numPermittedErrs),
                                        const int goodAlignmentScore,
                                        SegmentList& segments) const {
    // Verify consistency
    ASSERT ( refBatonList.getNmerSize() == othBatonList.getNmerSize() );
    // Compute all possible encoding values for the n-mers
    const int MaxNmerValue = (1 << (refBatonList.getNmerSize() * 2));
    // For each baton pair that falls within the respective windows in
    // the two baton lists, track the best alignment thus far..
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
                // refBaton and othBaton by forming an initial
                // segment.
                Segment initSeg =
                    makeSegment(refSeq, refBaton->getStartIndex(), 0,
                                refSeq.size(), otherSeq,
                                othBaton->getStartIndex(), 0, otherSeq.size());
                // See if we have a long-enough initial segement
                if (initSeg.getLength() < goodAlignmentScore) {
                    // Did not find any major alignment here.
                    // Onto the next pair.
                    continue;
                }
                // Found a major alignment. Build more sub-alignment
                // segments to the left and right of the initial
                // segment.
                int score = makeMultipleSegments(refSeq, codedRefNmers,
                                                 otherSeq, codedOthNmers,
                                                 initSeg, segments);
                return score;
            }
        }
    } // n-mer for-loop
    
    // When control drops here, we did not find a good-enough alignment
    return -1;
}

// This is a custom helper method that is used to search to the left
// and to the right of an identical baton-pair starting positions to
// locate first position where two consecutive nucleotide differences
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
    return Segment(baton1Pos - leftSpanLen, baton2Pos - leftSpanLen,
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
    ASSERT( baton1Pos < (int) refSeq.size() );
    ASSERT( baton2Pos < (int) otherSeq.size() );    
    // Ensure we can even move to the left if possible.
    if ((baton1Pos <= refSeqStartPos) || (baton2Pos <= othSeqStartPos)) {
        // The baton positions are already at or before start of
        // sequences
        return 0;
    }
    // Setup variables used in the do..while loop below to determine
    // number of different nucleotides seen thus far.
    register int refSeqPos = baton1Pos;  // index in refSeq
    register int othSeqPos = baton2Pos;  // index in otherSeq
    register int diffNtCnt = 0;          // Count of different NTs thus far

    // This loop typically ends when either currRefPos or currOthPos
    // becomes 0 or when diffNtCnt == 2.
    do {
        // Check the previous nucleotides for (mis)match
        refSeqPos--;
        othSeqPos--;
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
                refSeqPos += 2;
                othSeqPos += 2;
                break;
            }
        }
    } while ((refSeqPos > refSeqStartPos) && (othSeqPos > othSeqStartPos));
    
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
    register int othSeqPos = baton2Pos + 1;  // index int otherSeq
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
                refSeqPos--;
                othSeqPos--;
                break;
            }
        }
        // Check the next pair of nucleotides for (mis)match
        refSeqPos++;
        othSeqPos++;
    }
    // Return number of positions moved to right in reference sequence
    // (in case you are wondering, it will be the same in other
    // sequence as well). However, we never return a negative value.
    return refSeqPos - baton1Pos;
}

int
DefaultSequenceAligner::makeMultipleSegments(const std::string& refSeq,
                                             const IntVector& codedRefNmers,
                                             const std::string& otherSeq,
                                             const IntVector& codedOthNmers,
                                             const Segment& refSeg,
                                             SegmentList& segments) const {
    // Clear out any segments in the segments list
    segments.clear();
    // Reserve memory to minimize reallocs.
    segments.reserve(20);

    // Variable to track the total number of nucleotides aligned which
    // is sum of the segement length.
    int sumSegLen = 0;
    
    // First form segments to the left of the reference sequence.
    Segment currSeg = refSeg; 
    int refSeedPos  = -1;
    int othSeedPos  = -1;
    // Repeatedly form segments to the left of the current segement
    // using a seed n-mer.
    while (findSeedToLeft(codedRefNmers, currSeg.getRefESTOffset(),
                          codedOthNmers, currSeg.getOthESTOffset(),
                          refSeedPos, othSeedPos)) {
        // Form a new segment around the Seeds
        currSeg = makeSegment(refSeq, refSeedPos, 0,
                              currSeg.getRefESTOffset(),
                              otherSeq, othSeedPos, 0,
                              currSeg.getOthESTOffset());
        // Add segment to the vector (albeit in reverse order which we
        // will need to fix below).
        segments.push_back(currSeg);
        // Track number of nucleotides aligned
        sumSegLen += currSeg.getLength();
    }
    // Since in the above while-loop we added segments in a
    // right-to-left order we need to reverse them here to get a
    // left-to-right order.
    std::reverse(segments.begin(), segments.end());

    // Add the reference segment to the list now that we know all
    // segements that are to the left of the reference segment.
    segments.push_back(refSeg);
    sumSegLen += refSeg.getLength();
    
    // Repeatedly form segments to the right of the current segement
    // using a seed n-mer.
    currSeg = refSeg;
    while (findSeedToRight(codedRefNmers, currSeg.getRefESTEndOffset(),
                           codedOthNmers, currSeg.getOthESTEndOffset(),
                           refSeedPos, othSeedPos)) {
        // Form a new segment around the Seeds
        currSeg = makeSegment(refSeq, refSeedPos,
                              currSeg.getRefESTEndOffset() + 1,
                              refSeq.size(), otherSeq, othSeedPos,
                              currSeg.getOthESTEndOffset() + 1,
                              otherSeq.size());
        // Add segment to the vector (here it is in the correct order)
        segments.push_back(currSeg);
        // Track number of nucleotides aligned
        sumSegLen += currSeg.getLength();
    }   
    // Return total number of nucleotides aligned
    return sumSegLen;
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
                                       int& refSeedPos,
                                       int& othSeedPos) const {
    // Compute the left-position where we are going to stop searching
    const int refStopPos = std::max(0, refStartPos - SEED_SEARCH_DIST);
    const int othStopPos = std::max(0, othStartPos - SEED_SEARCH_DIST);
    // Check to see if we have something to search for
    if ((refStartPos <= 0) || (othStartPos <= 0)) {
        // No seeds to search for.
        return false;
    }
    // Find position of matching n-mers (aka seeds) from the two coded
    // vectors that are in nearest proximity to each other.
    return findNearestSeed(codedRefNmers, refStopPos, refStartPos - 1,
                           codedOthNmers, othStopPos, othStartPos - 1,
                           refSeedPos, othSeedPos);

    /*
    // Find position of matching n-mers (aka seeds) from the two coded
    // vectors.
    for(int refPos = refStartPos -1; (refPos >= refStopPos); refPos--) {
        const int refNmer = codedRefNmers[refPos];
        for(int othPos = othStartPos - 1; (othPos >= othStopPos); othPos--) {
            if (refNmer == codedOthNmers[othPos]) {
                // Found matching n-mers. We are done searching!
                refSeedPos = refPos;
                othSeedPos = othPos;
                return true;
            }
        }
    }
    // When control drops here that means we did not find a matching
    // seed.
    return false;
    */
}

// Helper method to determine start (or seed) for a segment to the
// right of a given pair (ref for reference sequence, and oth for
// other sequence being aligned with ref) of positions in two cDNA
// fragments.
bool
DefaultSequenceAligner::findSeedToRight(const IntVector& codedRefNmers,
                                        const int refStartPos,
                                        const IntVector& codedOthNmers,
                                        const int othStartPos,
                                        int& refSeedPos,
                                        int& othSeedPos) const {
    // Compute the left-position where we are going to stop searching
    const int refStopPos = std::min((int) codedRefNmers.size(),
                                    refStartPos + SEED_SEARCH_DIST);
    const int othStopPos = std::min((int) codedOthNmers.size(),
                                    othStartPos + SEED_SEARCH_DIST);
    // Check to see if we have something to search for
    if ((refStartPos >= refStopPos) || (othStartPos >= othStopPos)) {
        // No seeds to search for.
        return false;
    }
    // Find position of matching n-mers (aka seeds) from the two coded
    // vectors that are in nearest proximity to each other.
    return findNearestSeed(codedRefNmers, refStartPos + 1, refStopPos,
                           codedOthNmers, othStartPos + 1, othStopPos,
                           refSeedPos, othSeedPos);
}

bool
DefaultSequenceAligner::findNearestSeed(const IntVector& codedRefNmers,
                                        const int refStartPos,
                                        const int refEndPos,
                                        const IntVector& codedOthNmers,
                                        const int othStartPos,
                                        const int othEndPos,
                                        int& refSeedPos,
                                        int& othSeedPos) const {
    ASSERT ( refStartPos <= refEndPos );
    ASSERT ( othStartPos <= othEndPos );
    // Create a sorted array of SeedInfo objects for the reference
    // sequence.
    std::vector<SeedInfo> refSeedInfo;
    for(int ref = refStartPos; (ref < refEndPos); ref++) {
        refSeedInfo.push_back(std::make_pair(codedRefNmers[ref], ref));
    }
    // Sort the reference seed information for rapid lookups.
    std::sort(refSeedInfo.begin(), refSeedInfo.end(), LessSeedInfo());
    // Find the seed in the other sequence that has the smallest
    // distance to the reference sequence.
    bool matchFound  = false;
    int  minDistance = std::numeric_limits<int>::max(); // Variable used below
    for(int oth = othStartPos; ((oth < othEndPos) && (minDistance > 0)); oth++) {
        const int othNmer = codedOthNmers[oth];
        // Binary search in refSeedInfo for first matching value.
        std::vector<SeedInfo>::iterator matchIter =
            std::lower_bound(refSeedInfo.begin(), refSeedInfo.end(),
                             std::make_pair(othNmer, oth), LessSeedInfo());
        // For any matches found, track the closest one.
        while ((matchIter != refSeedInfo.end()) &&
               (matchIter->first == othNmer)) {
            // Track the closest matching n-mer thus far.
            const int distance = abs((matchIter->second - refStartPos) -
                                     (oth - othStartPos));
            if (minDistance > distance) {
                // Track the closest match.
                matchFound  = true;
                refSeedPos  = matchIter->second;
                othSeedPos  = oth;
                minDistance = distance;
            }
            // On to the next matching entry.
            matchIter++;
        }
    }
    // Return overall result as to whether matching seeds were found
    return matchFound;
}

#endif
