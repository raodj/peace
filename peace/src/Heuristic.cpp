#ifndef HEURISTIC_CPP
#define HEURISTIC_CPP

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

#include "Heuristic.h"
#include "Utilities.h"

Heuristic::Heuristic(const std::string& name)
    : refESTidx(-1), heuristicName(name), runCount(0), successCount(0) {
    // Nothing else to be done for now.
}

Heuristic::~Heuristic() {
    // Empty constructor begets an empty destructor
}

bool
Heuristic::shouldAnalyze(const int otherEST) {
    bool value;
    if ((value = runHeuristic(otherEST)) == true) {
        // This pair should be analyzed further. 
        successCount++;
    }
    // Track number of times this heuristic was run.
    runCount++;
    // Return true (indicating further analysis is needed) or false
    // (no further analysis required).
    return value;
}

void
Heuristic::printStats(std::ostream& os) const {
    os << "Statistics from " << getName()   << " heuristic:\n"
       << "\tNumber of calls        : " << runCount     << "\n"
       << "\tNumber of successes    : " << successCount << std::endl;
}

#endif
