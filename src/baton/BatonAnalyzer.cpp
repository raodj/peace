#ifndef BATON_ANALYZER_CPP
#define BATON_ANALYZER_CPP

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

#include "BatonAnalyzer.h"
#include "AlignmentInfo.h"
#include <limits>

// The static variable for threshold
int BatonAnalyzer::windowSize = 100;
// The static variable for baton head size
int BatonAnalyzer::nMerSize   = 3;
// The static variable for threshold
int BatonAnalyzer::threshold  = 7;
// Acceptable number of differences when exploring left and right
int BatonAnalyzer::numPermittedErrs = 10;
// The good-enough alignment score command line argument
int BatonAnalyzer::goodAlignmentScore = 15;

// The set of command line arguments used by this analyzer.
arg_parser::arg_record BatonAnalyzer::argsList[] = {
    {"--threshold", "Min. baton count above which 2 reads are similar (def=7)",
     &BatonAnalyzer::threshold, arg_parser::INTEGER},
    {"--nMers", "The number of bases (1 to 3) to be used for baton ends",
     &BatonAnalyzer::nMerSize, arg_parser::INTEGER},
    {"--window", "The average size (in nt) of a window for matching identical batons",
     &BatonAnalyzer::windowSize, arg_parser::INTEGER},
    {"--permErrs", "#of differences to be accepted when checking candidate alignments",
     &BatonAnalyzer::numPermittedErrs, arg_parser::INTEGER},        
    {"--goodScore", "A good-enough alignment score to short circuit exhaustive analysis",
     &BatonAnalyzer::goodAlignmentScore, arg_parser::INTEGER},    
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

BatonAnalyzer::BatonAnalyzer(const int refESTidx,
                             const std::string& outputFileName)
    : ESTAnalyzer("baton", refESTidx, outputFileName) {
    // Clear out the reference sequence baton information
    refBatonList = NULL;
}

BatonAnalyzer::~BatonAnalyzer() {
    // Free up the memory for the baton lists maintained in the cache.
    
    for(size_t idx = 0; (idx < blCache.size()); idx++) {
        if (blCache[idx] != NULL) {
            delete blCache[idx];
        }
    }
}

// Interface/API method defined in ESTAnalyzer base class
void
BatonAnalyzer::showArguments(std::ostream& os) {
    ESTAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(BatonAnalyzer::argsList);
    os << ap;
}

// Interface/API method defined in ESTAnalyzer base class
bool
BatonAnalyzer::parseArguments(int& argc, char **argv) {
    // Let the overloaded method do the actual task.
    return parseArguments(argc, argv, false);
}

bool
BatonAnalyzer::parseArguments(int& argc, char **argv,
                              const bool assemblerFlag) {
    arg_parser ap(BatonAnalyzer::argsList);
    ap.check_args(argc, argv, false);
    // Ensure threshold is valid.
    if (threshold < 1) {
        std::cerr << analyzerName
                  << ": Threshold must be >= 1 (suggested value: 7)"
                  << "(use --threshold option)\n";
        return false;
    }
    // Check and update the goodalignmentScore
    if (goodAlignmentScore == -1) {
        goodAlignmentScore = std::numeric_limits<int>::max();
    }
    if (assemblerFlag) {
        // In this case, don't call base class methods as this class
        // is being used for assembly.
        return true;
    }
    // Now let the base class do processing and return the result.
    return ESTAnalyzer::parseArguments(argc, argv);
}

// Interface/API method defined in ESTAnalyzer base class
int
BatonAnalyzer::initialize() {
    // Use the overloaded method to get the job done.
    return initialize(true);
}

int
BatonAnalyzer::initialize(const bool loadData) {
    if (loadData) {
        // First load the data based on the file specified.
        if ((estFileName != NULL) && (!loadFASTAFile(estFileName))) {
            // Loading EST's from FASTA file was not successful.  Can't do
            // much further.
            return 1;
        }
        if ((sffFileName != NULL) && (!loadSFFFile(sffFileName))) {
            // Loading EST's from SFF file was not successful.  Can't do
            // much further.
            return 2;
        }
    }
    // Now that we know the number of cDNA fragments to be processed,
    // initialize our cache with NULL values. We reserve twice the
    // number of entries for normal and reverse complement baton list
    blCache.reserve(EST::getESTCount() * 2);
    blCache.resize (EST::getESTCount() * 2, NULL);
    
    // Check and warn that heuristics are not really being used, if we
    // have a heuristic chain specified.
    if (chain != NULL) {
        std::cerr << analyzerName
                  << ": This analyzer does not use any heuristics. So "
                  << "the heuristic chain is being ignored.\n";
    }
    // Initialization successful.
    return 0;
}

// Interface/API method defined in ESTAnalyzer base class
int
BatonAnalyzer::setReferenceEST(const int estIdx) {
    // Clear out any baton list created for a direct-given consensus
    // type sequence.
    if ((refESTidx == -1) && (refBatonList != NULL)) {
        delete refBatonList;
        refBatonList = NULL;
    }
    
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        // Pre-compute normal baton list, if it is not already in our
        // blCache...
        refBatonList = getBatonList(estIdx, false);
        ASSERT ( refBatonList != NULL );
        // Build up the n-mer list for future use.
        refBatonList->getCodedNMers(codedRefNmers);
        // Save reference EST idx for future use.
        refESTidx = estIdx;
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

int
BatonAnalyzer::setReferenceEST(const char* refSeq) {
    // Clear out any earlier baton list created for a direct-given
    // consensus type sequence.
    if ((refESTidx == -1) && (refBatonList != NULL)) {
        delete refBatonList;
        refBatonList = NULL;
    }
    // Pre-compute normal baton list for the given reference sequence.
    refBatonList = new BatonList(refSeq, nMerSize, false, windowSize);
    ASSERT ( refBatonList != NULL );
    // Build up the n-mer list for future use.
    refBatonList->getCodedNMers(codedRefNmers);
    // Set reference EST idx to an invalid value for consistency and
    // cross referencing in the future.
    refESTidx = -1;
    return 0; // everything went well
}

// Implement interface method from ESTAnalyzer base class.
float
BatonAnalyzer::getMetric(const int otherEST) {
    // First try comparing the "normal" baton lists for the reference
    // and other cDNA fragment.
    const BatonList* const refBatonList = getBatonList(refESTidx, false);
    const BatonList* const othBatonList = getBatonList(otherEST,  false);
    ASSERT ( refBatonList != NULL );
    ASSERT ( othBatonList != NULL );
    // Compute the similarity metric.
    const float normMetric = getMetric(refBatonList, othBatonList);
    // Now obtain the reverse complement (RC) baton list for other
    const BatonList* const rcBatonList = getBatonList(otherEST,  true);
    // Comptue the RC similarity metric.
    const float rcMetric = getMetric(refBatonList, rcBatonList);
    // Return the best of the two as the final similarity metric
    // std::cout << "normMetric = " << normMetric
    //           << ", rcMetric = " << rcMetric << std::endl;
    return std::max(normMetric, rcMetric);
}

// Helper for the getMetric() method defined earlier.
float
BatonAnalyzer::getMetric(const BatonList* const list1,
                         const BatonList* const list2) const {
    ASSERT ( list1 != NULL );
    ASSERT ( list2 != NULL );
    // Get the high frequency window pairs..
    std::vector<WindowPair> pairList;
    list1->getWindowPairs(*list2, pairList, threshold);
    // Use the size of the pair list as the similarity
    return pairList.size();
}

int
BatonAnalyzer::getBestAlignmentInWindow(const BatonList& refBatonList,
                                        const IntVector& codedRefNmers,
                                        const int refWindow,
                                        int&  refAlignPos,
                                        const BatonList& othBatonList,
                                        const IntVector& codedOthNmers,
                                        const int othWindow,
                                        int&  othAlignPos,
                                        const int numPermittedErrs,
                                        const int goodAlignmentScore) const {
    // Verify consistency
    ASSERT ( refBatonList.getNmerSize() == nMerSize );
    ASSERT ( othBatonList.getNmerSize() == nMerSize );
    // Compute all possible encoding values for the n-mers
    const int MaxNmerValue = (1 << (nMerSize * 2));
    // For each baton pair that falls within the respective windows in
    // the two baton lists, track the best alignment thus far..
    int bestScore = 0;
    for(int nMer = 0; (nMer < MaxNmerValue); nMer++) {
        for(NmerBatonList::const_iterator refBaton = refBatonList[nMer].begin();
            (refBaton != refBatonList[nMer].end()); refBaton++) {
            if (!refBatonList.isBatonInWindow(*refBaton, refWindow)) {
                // This baton is not in the window are are
                // checking. So skip it.
                continue;
            }
            for(NmerBatonList::const_iterator othBaton = othBatonList[nMer].begin();
                (othBaton != othBatonList[nMer].end()); othBaton++) {
                if (!othBatonList.isBatonInWindow(*othBaton, othWindow)) {
                    // This baton is not in the window are
                    // checking. So skip it.
                    continue;
                }
                // Determine alignment score for the current pair of
                // refBaton and othBaton
                int currScore =
                    getAlignmentScore(codedRefNmers, refBaton->getStartIndex(),
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

int
BatonAnalyzer::getAlignmentScore(const IntVector& codedRefNmers, int baton1Pos,
                                 const IntVector& codedOthNmers, int baton2Pos,
                                 int numPermittedErrs) const {
    // Set the maximum length of alignment we can achieve based on the
    // shorter of the two cDNA fragments.
    const int MaxAlignLen = std::min(codedRefNmers.size(), codedOthNmers.size());
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
            tryMoveToLeft  = true;      // Try going left more in next iteration
            tryMoveToRight = false;     // Don't go right in next iteration
        }
        // When control drops here we have encountered one difference
        // either to the left or to the right. So track it.
        numPermittedErrs--;
    } while ((numPermittedErrs > 0) && ((leftMoves + rightMoves)< MaxAlignLen));
    // Return the sum of the left and right distances we have managed
    // to move while accruing numPermittedErrs number of differences/errors.
    return (leftMoves + rightMoves);
}

// Helper method to manage cached batons. All the methods that need
// access to a baton list should obtain it via this method.
BatonList*
BatonAnalyzer::getBatonList(const int estIdx, const bool getRC) {
    // Determine index position for the cached entry.  The normal and
    // RC baton lists for an cDNA fragment with index k are stored in
    // at position k*2 and k*2+1 respectively.
    const int cachePos = estIdx * 2 + (getRC ? 1 : 0);
    // If we don't have an entry at the corresponding position, then
    // we need to construct one now.
    if (blCache[cachePos] == NULL) {
        // Don't have a baton list here. So create one and save it for
        // further use.
        blCache[cachePos] = new BatonList(estIdx, nMerSize, getRC, windowSize);
        // Print the baton list we just built for testing...
        // std::cout << "EST #" << estIdx << std::endl;
        // std::cout << *blCache[cachePos] << std::endl;
    }
    // When control drops here we always have a valid entry in the cache
    return blCache[cachePos];
}

// Interface method from ESTAnalyzer
bool
BatonAnalyzer::getAlignmentData(int &alignmentData) {
    alignmentData = -1;
    return false;
}

// Primary method method to compute the alignment between the reference
// sequence and a given EST.
bool
BatonAnalyzer::align(const int otherEST, AlignmentInfo& info) {
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
        info = AlignmentInfo(otherEST, refESTidx, refAlignPos - othAlignPos,
                             score, false);
        return true;
    }
    // Since first round of analysis failed, next check to see if
    // alignment with reverse complement) nucleotide sequence of
    // otherEST yeilds good results.
    if (align(otherEST, refAlignPos, othAlignPos, score, true)) {
        // Found a sufficiently good alignment. Update info and return
        info = AlignmentInfo(otherEST, refESTidx, refAlignPos - othAlignPos,
                             score, true);
        return true;
    }
    // Both normal and reverse complement analysis failed.
    return false;
}

// Helper method method to compute the alignment between the reference
// sequence and a given EST.
bool
BatonAnalyzer::align(const int otherEST, int& refAlignPos, int& othAlignPos,
                     int& score, const bool doRC) {
    // We should have a reference baton list already set
    ASSERT( refBatonList != NULL );
    // Build the normal (not reverse complement) baton list for the
    // otherEST for analysis.
    const BatonList *othBatonList = getBatonList(otherEST, doRC);
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
    // Set starting default score to zero
    score = 0;
    // Search among pairs of windows returned in the list for best
    // alignment.  But first get the coded n-mers for the otherEST.
    IntVector codedOthNmers;
    othBatonList->getCodedNMers(codedOthNmers);
    for(size_t winIdx = 0; (winIdx < pairList.size()); winIdx++) {
        // Get alignment based on the current window-pair
        int currRefAlignPos, currOthAlignPos;
        int currScore =
            getBestAlignmentInWindow(*refBatonList, codedRefNmers,
                                     pairList[winIdx].first, currRefAlignPos,
                                     *othBatonList, codedOthNmers,
                                     pairList[winIdx].second, currOthAlignPos,
                                     numPermittedErrs, goodAlignmentScore);
        // Track the best alignment we have had so far...
        if (currScore > score) {
            score       = currScore;
            refAlignPos = currRefAlignPos;
            othAlignPos = currOthAlignPos;
        }
        // Short circuit exhaustive alignment search if a reasonably
        // good alignment was found.
        if (score >= goodAlignmentScore) {
            break;
        }
    }
    // Indicate we have an alignment
    return true;
}

int
BatonAnalyzer::analyze() {
    // This method needs to be implemented in a fashion that is
    // similar to FWAnalyzer::analyze() method.
    return 0;
}

#endif
