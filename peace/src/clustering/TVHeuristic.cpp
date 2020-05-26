#ifndef TV_HEURISTIC_CPP
#define TV_HEURISTIC_CPP

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

#include "TVHeuristic.h"
#include "ESTCodec.h"
#include "ESTList.h"
#include "ArgParser.h"
#include "HeuristicChain.h"
#include "RuntimeContext.h"
#include "ParameterSetManager.h"

TVHeuristic::TVHeuristic(HeuristicChain *chain) : NewUVHeuristic(chain, "tv") {
    uvSuccessCount = 0;
    t              = 65;
    windowLen      = 100;
    // Set up current adaptive parameter set index to an invalid value.
    currParamSetIndex = -1;    
}

TVHeuristic::~TVHeuristic() {
    // matchTable is now a vector and no longer needs to be cleaned-up
    // here.
}

void
TVHeuristic::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add its parameters first
    NewUVHeuristic::addCommandLineArguments(argParser);
    // The set of arguments for this class.  Note that some of the base
    // class static instance variables are reused here so that the values
    // are consistently set.
    const ArgParser::ArgRecord ArgsList[] = {
        {"--tv_t", "t (number of v-word matches)",
         &t, ArgParser::INTEGER},
        {"--tv_win", "Window size for t/v heuristics",
         &windowLen, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };    
    // Add additional command line arguments for the user to work with.
    argParser.addValidArguments(ArgsList);
}

bool
TVHeuristic::initialize() {
    // Validate our command-line parameter(s).
    if (t < 1) {
        std::cerr << heuristicName
                  << ": t (number of v-words matches) >= 1 "
                  << "(use --tv_t option)\n";
        return false;
    }
    
    // First let base class do initialization
    if (!NewUVHeuristic::initialize()) {
        // Error initializing base class
        return false;
    }
    
    // First compute the longest EST we have.  For this we need to get
    // hold of the ESTList from the RuntimeContext object in the
    // SubSystem.
    ASSERT ( heuristicChain != NULL );
    ASSERT ( heuristicChain->getESTList() != NULL );
    const ESTList* estList = heuristicChain->getESTList();
    size_t maxESTlen = estList->getMaxESTLen();
    int maxFrameSize = heuristicChain->getParamSetMgr()->getMaxFrameSize();
    // Add extra windowLen characters to ease processing.
    matchTable.resize(maxESTlen + maxFrameSize + v, 0);
    // Everything went well
    return true;
}

int
TVHeuristic::setReferenceEST(const EST* est) {
    ASSERT( est != NULL );
    if ((refEST != NULL) && (refEST->getID() == est->getID())) {
        // The current reference EST is the same as parameter. Nothing
        // to be done.
        return 0;
    }
    // Simply let the base class to its job
    return NewUVHeuristic::setReferenceEST(est);
}

bool
TVHeuristic::runHeuristic(const EST* otherEST) {
    ASSERT( otherEST != NULL );
    // Call updateParameters method to use proper heuristic params
    if (!updateParameters(otherEST)) {
        // This pair need not be analyzed further.
        return false;
    }
    // First check if this pair passes UV Sample heuristic.  If not do
    // no further analysis.
    if (!NewUVHeuristic::runHeuristic(otherEST)) {
        // This pair need not be analyzed further.
        return false;
    }    
    // Track number of successful base-class checks.
    uvSuccessCount++;
    // Now apply tv-heuristic to see if this pair should be analyzed
    // further.
                     
    int numMatches = 0;
    // bitsToShift is set to 2*(v-1) in the base class.
    ESTCodec::NormalEncoder<bitsToShift, BitMask> encoder;
    if (bestMatchIsRC) {
        numMatches = countCommonWords(otherEST, encoder, s1RCWordMap);
    } else {
        numMatches = countCommonWords(otherEST, encoder, s1WordMap);
    }
    // Ensure number of matches exceeds threshold limits
    return (numMatches >= t);
}

bool
TVHeuristic::updateParameters(const EST* otherEST) {
    const ParameterSetManager* const paramMgr =
        heuristicChain->getParamSetMgr();
    ASSERT( paramMgr != NULL );
    const int paramSetIndex =
        paramMgr->getParameterSet(refEST->getID(), otherEST->getID());
    if (paramSetIndex == -1) {
        // These two fragments are very different in lengths. Don't
        // bother comparing them at all.
        return false;
    }
    if (paramSetIndex == currParamSetIndex) {
        // The parameters that we are currently using are just
        // fine. No need to update them any further. Analyze the two
        // fragments using the current parameter set.
        return true;
    }
    // Update and move to a new set of parameters.
    const ParameterSet* const parameterSet =
        paramMgr->getParameterSet(paramSetIndex);
    // Assign parameters according to parameter set chosen
    windowLen = parameterSet->frameSize;
    t         = parameterSet->t;
    u         = parameterSet->u;
    wordShift = parameterSet->wordShift;
    passes    = (u < 6) ? 2 : 3;
    // Save the current parameter index for future short circuiting
    currParamSetIndex = paramSetIndex;
    // Analyze the two fragments using the current parameter set.
    return true;
}

void
TVHeuristic::printStats(std::ostream& os) const {
    // First let the base class do its thing.
    NewUVHeuristic::printStats(os);
    // Display additional information about uv success
    os << "\tNumber of u/v successes: " << uvSuccessCount << std::endl;
}

#endif
