#ifndef HEURISTIC_CPP
#define HEURISTIC_CPP

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

#include "Heuristic.h"
#include "Utilities.h"

Heuristic::Heuristic(const std::string& name, HeuristicChain* chain)
    : Component(name), refEST(NULL), heuristicName(name),
      heuristicChain(chain), runCount(0),
      successCount(0) {
    ASSERT( chain != NULL );
    // Nothing else to be done for now.
}

Heuristic::~Heuristic() {
    // Empty constructor begets an empty destructor
}

bool
Heuristic::shouldAnalyze(const EST* otherEST) {
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
