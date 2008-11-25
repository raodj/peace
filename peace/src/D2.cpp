#ifndef D2_CPP
#define D2_CPP

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

#include "D2.h"
#include "EST.h"

D2::D2(const int refESTidx, const std::string& outputFileName)
    : ESTAnalyzer("D2", refESTidx, outputFileName) {
}

D2::~D2() {
}

void
D2::showArguments(std::ostream& os) {
    ESTAnalyzer::showArguments(os);
    // Currently this class does not have any specific command line
    // arguments.  Consequently, there is nothing further to be done.
}

bool
D2::parseArguments(int& argc, char **argv) {
    // Currently this class does not have any specific command line
    // arguments.  The base class does all the necessary tasks.

    // Let the base class do processing and return the result.
    return ESTAnalyzer::parseArguments(argc, argv);
}

int
D2::initialize() {
    // Currently this class does not require any special
    // initialization to be performed.
    return 0;
}

int
D2::setReferenceEST(const int estIdx) {
    if ((estIdx >= 0) && (estIdx < EST::getESTCount())) {
        refESTidx = estIdx;
        return 0;
    }
    // Invalid est index.
    return 1;
}

float
D2::analyze(const int otherEST) {
    if ((otherEST >= 0) && (otherEST < EST::getESTCount())) {
        return 0;
    }
    // Invalid est index.
    return -1;
}

int
D2::analyze() {
    return -1;
}

#endif
