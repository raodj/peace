#ifndef EST_ANALYZER_CPP
#define EST_ANALYZER_CPP

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

#include "ESTAnalyzer.h"
#include "MPIHelper.h"
#include "ArgParser.h"
#include "SubSystem.h"
#include "RuntimeContext.h"
#include "HeuristicChain.h"
#include "EST.h"

ESTAnalyzer::ESTAnalyzer(const std::string& name) :
    Component(name), refEST(NULL), chain(NULL), estList(NULL),
    initialRefESTidx(0) {
    //  Nothing else to be done.
}

ESTAnalyzer::~ESTAnalyzer() {
    // Empty constructor begets an empty destructor
}

int
ESTAnalyzer::setHeuristicChain(HeuristicChain* hChain) {
    chain = hChain;
    return 0; // Everything went well
}

void
ESTAnalyzer::setESTList(ESTList* estList) {
    this->estList = estList;
}

bool
ESTAnalyzer::initialize() {
    // Do base class initialization first
    if (!Component::initialize()) {
        // Base class initialization failed!
        return false;
    }
    // Ensure the convenience pointer to the shared list of ESTs is
    // set.
    if (estList == NULL) {
        std::cerr << "Error: ESTAnalyzer does not have a valid list of ESTs "
                  << "to use.\n"
                  << "Ensure ESTAnalyzer::setESTList() has been called.\n";
        return false;
    }
    // Check and set EST list for heuristic chains as needed
    if ((chain != NULL) && (chain->getESTList() == NULL)) {
        // Setup the same EST list with the heuristic chain.
        chain->setESTList(estList);
    }
    // Ensure the heuristics and the analyzer are operating with the
    // same set of data.
    if ((chain != NULL) && (chain->getESTList() != estList)) {
        std::cerr << "Error: It appears that the ESTAnalyzer "
                  << getName() << " and the heuristics (associated with "
                  << "the analyzer) are not operating on the same ESTList.\n";
        return false;
    }
    // Initialize the heuristic chain preparing it for use during
    // clustering.
    if ((chain != NULL) && (!chain->initialize())) {
        // Error occured during initialization. Bail out.
        return false;
    }
    // Initialization successful.
    return true;
}

void
ESTAnalyzer::finalize() {
    estList = NULL;
}

float
ESTAnalyzer::analyze(const EST* otherEST, const bool useHeuristics,
                     const bool useHeavyWeight) {
    // Check with the heuristic chain
    if (useHeuristics && (chain != NULL) &&
        !chain->shouldAnalyze(otherEST)) {
        // Heuristics indicate we should not do D2. So skip it.
        return getInvalidMetric();
    }
    if (useHeavyWeight) {
        return getMetric(otherEST);
    }
    // Indicate that we would have used heavy weight analysis
    return getValidMetric();
}

void
ESTAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    const ArgParser::ArgRecord CommonArgsList[] = {
        {"--extIdx", "Index of the initial reference EST",
         &initialRefESTidx, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    argParser.addValidArguments(CommonArgsList);
}

ESTAnalyzer& 
ESTAnalyzer::operator=(const ESTAnalyzer&) {
    return *this;
}

void
ESTAnalyzer::displayStats(std::ostream& os) {
    if (chain != NULL) {
        chain->printStats(os, MPI_GET_RANK());
    }
}

#endif
