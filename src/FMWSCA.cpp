#ifndef FMWSCA_CPP
#define FMWSCA_CPP

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

#include "FMWSCA.h"

// The static instance variables for command line arguments.
bool FMWSCA::caseSensitive  = false;

// The common set of arguments for all FW EST analyzers
arg_parser::arg_record FMWSCA::argsList[] = {
    {"--case", "Make comparisons of BP case sensitive",
     &FMWSCA::caseSensitive, arg_parser::BOOLEAN},
    {NULL, NULL, NULL, arg_parser::BOOLEAN}
};

FMWSCA::FMWSCA(const int refESTidx, const std::string& outputFile)
    : FWAnalyzer(std::string("fmwsca"), refESTidx, outputFile) {
    // Nothing else to be done for now.
}

FMWSCA::~FMWSCA() {
    // Nothing else to be done for now.
}

void
FMWSCA::showArguments(std::ostream& os) {
    FWAnalyzer::showArguments(os);
    // Use a arg parser object to conveniently display common options.
    arg_parser ap(FMWSCA::argsList);
    os << ap;
}

bool
FMWSCA::parseArguments(int& argc, char **argv) {
    arg_parser ap(FMWSCA::argsList);
    ap.check_args(argc, argv, false);
    // Now let the base class do processing and return the result.
    return FWAnalyzer::parseArguments(argc, argv);
}

float
FMWSCA::getMetric(const std::string& refFrame,
                 const std::string& otherFrame, const int wordSize) {
    const int LastWordPos = (int) refFrame.length() - wordSize;
    int matchCount        = 0;
    for(int wordPos = 0; (wordPos < LastWordPos); wordPos++) {
        const std::string word = refFrame.substr(wordPos, wordSize);
        int searchPos = 0;
        do {
            if ((searchPos = (int) otherFrame.find(word, searchPos)) != -1) {
                matchCount++;
                searchPos++;
            }
        } while (searchPos != -1);
    }

    return ((float) matchCount) / (LastWordPos * LastWordPos);
}

#endif
