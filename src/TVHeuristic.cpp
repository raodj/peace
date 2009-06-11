#ifndef TV_HEURISTIC_CPP
#define TV_HEURISTIC_CPP

//---------------------------------------------------------------------------
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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "TVHeuristic.h"
#include "ESTCodec.h"
#include "EST.h"

// Define the static parameters
int TVHeuristic::t = 65;

int TVHeuristic::windowLen = 100;

// The set of arguments for this class.  Note that some of the base
// class static instance variables are reused here so that the values
// are consistently set.
arg_parser::arg_record TVHeuristic::argsList[] = {
    {"--tv_t", "t (number of v-word matches) (default=5)",
     &TVHeuristic::t, arg_parser::INTEGER},
    {NULL, NULL}
};

TVHeuristic::TVHeuristic(const std::string& outputFileName)
    : NewUVHeuristic("tv", outputFileName) {
    matchTable     = NULL;
    uvSuccessCount = 0;
}

TVHeuristic::~TVHeuristic() {
    if (matchTable != NULL) {
        delete [] matchTable;
    }
}

void
TVHeuristic::showArguments(std::ostream& os) {
    // First let the base class do the operation
    NewUVHeuristic::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(TVHeuristic::argsList);
    os << ap;
}

bool
TVHeuristic::parseArguments(int& argc, char **argv) {
    // Let the base class do the same operation first
    if (!NewUVHeuristic::parseArguments(argc, argv)) {
        // Error parsing base class parameters
        return false;
    }
    // Let's process parameters now.
    arg_parser ap(TVHeuristic::argsList);
    ap.check_args(argc, argv, false);
    // Ensure values are valid.
    if (t < 1) {
        std::cerr << heuristicName
                  << ": t (number of v-words matches) >= 1 "
                  << "(use --tv_t option)\n";
        return false;
    }
    // Everything went well
    return true;
}

int
TVHeuristic::initialize() {
    // First let base class do initialization
    NewUVHeuristic::initialize();
    // First compute the longest EST we have.
    size_t maxESTlen = EST::getMaxESTLen();
    
    // Add extra 100 characters to ease processing.
    matchTable = new char[maxESTlen + windowLen];
    // Everything went well
    return 0;
}

int
TVHeuristic::setReferenceEST(const int estIdx) {
    if (refESTidx == estIdx) {
        // The reference EST is the same. Nothing to be done.
        return 0;
    }
    // Simply let the base class to its job
    return NewUVHeuristic::setReferenceEST(estIdx);
}

bool
TVHeuristic::runHeuristic(const int otherEST) {
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
    // Print intermediate stats to compare with wcd
    // std::cout << "tv" << (bestMatchIsRC ? "'" : "")
    //           << "("  << refESTidx << ", " << otherEST << ") = "
    //           << numMatches << std::endl;
    
    // Ensure number of matches exceeds threshold limits
    return (numMatches >= TVHeuristic::t);
}

void
TVHeuristic::printStats(std::ostream& os) const {
    // First let the base class do its thing.
    NewUVHeuristic::printStats(os);
    // Display additional information about uv success
    os << "\tNumber of u/v successes: " << uvSuccessCount << std::endl;
}

#endif
