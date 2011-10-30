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
#include "BatonAlignmentInfo.h"
#include "ArgParser.h"
#include <limits>

BatonAnalyzer::BatonAnalyzer(const bool usedForAssembly) 
    : ESTAnalyzer("baton"), usedByAssembler(usedForAssembly) {
    // Initialize command-line arguments to default values
    threshold          = 7;
    numPermittedErrs   = 10;
    goodAlignmentScore = 15;
    // Clear out the reference sequence baton information
    refBatonList = NULL;
}

BatonAnalyzer::~BatonAnalyzer() {
}

void
BatonAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    ESTAnalyzer::addCommandLineArguments(argParser);
    // Now add our custom parameters.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--threshold", "Min. baton count above which 2 reads are similar",
         &threshold, ArgParser::INTEGER},
        {"--permErrs", "#of differences to be accepted when checking candidate alignments",
         &numPermittedErrs, ArgParser::INTEGER},        
        {"--goodScore", "A good-enough alignment score to short circuit exhaustive analysis",
         &goodAlignmentScore, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(ArgsList);
}

// Interface/API method defined in ESTAnalyzer base class
bool
BatonAnalyzer::initialize() {
    // Let the base class do its initialization
    if (!ESTAnalyzer::initialize()) {
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
    // Check and warn that heuristics are not really being used, if we
    // have a heuristic chain specified.
    if (chain != NULL) {
        std::cerr << getName()
                  << ": This analyzer does not use any heuristics. So "
                  << "the heuristic chain is being ignored.\n";
    }
    // Initialization successful.
    return true;
}

// Interface/API method defined in ESTAnalyzer base class
int
BatonAnalyzer::setReferenceEST(const EST* est) {
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
        refBatonList = blCache.getBatonList(est, false);
        ASSERT ( refBatonList != NULL );
        // Build up the n-mer list for future use.
        refBatonList->getCodedNMers(codedRefNmers);
        return 0; // everything went well
    }
    // Invalid est index.
    return 1;
}

int
BatonAnalyzer::setReferenceEST(const char* refSeq) {
    // Clear out any earlier baton list created for a direct-given
    // consensus type sequence.
    if (refBatonList != NULL) {
        delete refBatonList;
        refBatonList = NULL;
    }
    // Pre-compute normal baton list for the given reference sequence.
    refBatonList = new BatonList(refSeq, blCache.getBatonHeadSize(), false,
                                 blCache.getWindowSize());
    ASSERT ( refBatonList != NULL );
    // Build up the n-mer list for future use.
    refBatonList->getCodedNMers(codedRefNmers);
    // Set reference EST to an invalid value for consistency and cross
    // referencing in the future.
    refEST = NULL;
    return 0; // everything went well
}

// Implement interface method from ESTAnalyzer base class.
float
BatonAnalyzer::getMetric(const EST* otherEST) {
    ASSERT ( otherEST != NULL );
    ASSERT ( refEST   != NULL );
    // First try comparing the "normal" baton lists for the reference
    // and other cDNA fragment.
    const BatonList* const refBatonList = blCache.getBatonList(refEST, false);
    const BatonList* const othBatonList = blCache.getBatonList(otherEST,false);
    ASSERT ( refBatonList != NULL );
    ASSERT ( othBatonList != NULL );
    // Compute the similarity metric.
    const float normMetric = getMetric(refBatonList, othBatonList);
    // Now obtain the reverse complement (RC) baton list for other
    const BatonList* const rcBatonList = blCache.getBatonList(otherEST,  true);
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
    std::priority_queue<WindowPair> pairList;
    list1->getWindowPairs(*list2, pairList, threshold);
    // Use the highest frequency window count as a similarity metric
    return (float) (pairList.empty() ? 0 : pairList.top().getScore());
}

#endif
