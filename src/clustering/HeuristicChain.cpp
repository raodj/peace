#ifndef HEURISTIC_CHAIN_CPP
#define HEURISTIC_CHAIN_CPP

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

#include "HeuristicChain.h"
#include "HeuristicFactory.h"
#include "RuntimeContext.h"
#include "SubSystem.h"
#include "Utilities.h"

HeuristicChain::HeuristicChain() : Component("HeuristicChain"), estList(NULL) {
    // Set the last hint to an invalid value
    hints[LAST_HINT] = -1;
}

HeuristicChain::~HeuristicChain() {
    for(size_t i = 0; (i < chain.size()); i++) {
        delete chain[i];
    }
    chain.clear();
}

bool
HeuristicChain::addHeuristic(Heuristic* h) {
    ASSERT( h != NULL );
    chain.push_back(h);
    // Everything went well
    return true;
}

bool
HeuristicChain::initialize() {
    // Avoid redundant/duplicate initialization
    if (isInitialized()) {
        // Already initialized.
        return true;
    }
    // OK, first time initialization. First let the base class initialize
    if (!Component::initialize()) {
        // Base class initialization baild. Big boo boo
        return false;
    }
    // Inititalize the parameter set manager after setting up the
    // convenience reference to the shared estList for the parameter
    // set manager to use.
    if (estList == NULL) {
        std::cerr << "Error: HeuristicChain does not have a valid list of "
                  << "ESTs to use.\n"
                  << "Ensure HeuristicChain::setESTList() has been called.\n";
        return false;
    }
    psetMgr.setESTList(estList);
    psetMgr.initialize();
    // Setup reference to the parameter set manager to receive
    // notifications whenever the estlist changes
    if (subSystem != NULL) {
        ASSERT( subSystem->getContext() != NULL );
        subSystem->getContext()->setESTListListener(&psetMgr);
    }
    
    // Now inititalize each heuristic in the chain.
    for (size_t i = 0; (i < chain.size()); i++) {
	// Setup sub-system for each heuristic (in case it needs it)
	chain[i]->setSubSystem(subSystem);
        if (!chain[i]->initialize()) {
            // Error occured during initialization. Bail out.
            std::cerr << "Error initializing heuristic "
                      << chain[i]->getName() << std::endl;
            return false;
        }
    }
    return true; // Everything went well
}

int
HeuristicChain::setReferenceEST(const EST* est) {
    ASSERT( est != NULL );
    int result;
    for (size_t i = 0; (i < chain.size()); i++) {
        if ((result = chain[i]->setReferenceEST(est)) != 0) {
            // Error occured during setting up reference fragment. Bail out.
            return result;
        }
    }
    return 0; // Everything went well
}

Heuristic*
HeuristicChain::getHeuristic(const std::string& name) const {
    for(size_t i = 0; (i < chain.size()); i++) {
        if (chain[i]->getName() == name) {
            // Found a match!
            return chain[i];
        }
    }
    // A matching heuristic was not found.
    return NULL;
}

void
HeuristicChain::printStats(std::ostream& os, const int rank) const {
    os << "Heuristics on process with Rank " << rank << "\n"
       << "-------------------------------------------\n";
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->printStats(os);
    }
}


bool
HeuristicChain::setupChain(const std::string& heuristicStr) {
    // Check if there is anything to be created.
    if (heuristicStr.empty()) {
        return true;
    }
    // Process the non-empty heuristic string.
    std::string hStr(heuristicStr);
    // Process one word at a time. Words are separated by a "hyphen"
    // character.
    while (!hStr.empty()) {
        // Locate the next hyphen character and get name of heuristic
        const std::string::size_type hyphenLoc = hStr.find('-');
        const std::string name = hStr.substr(0, hyphenLoc);
        // Create heuristic and add it to the chain
        Heuristic *heuristic =
            HeuristicFactory::create(name, this);
        if (heuristic == NULL) {
            // Break out and return false, which will cause main to show usage
            return false;
        }
        // This assert caused errors if invalid heuristics were entered
        //ASSERT( heuristic != NULL );
        addHeuristic(heuristic);
        // Remove already created heuristic to process next one in chain.
        if (hyphenLoc == std::string::npos) {
            // No hyphen, so this was the last heuristic in the chain
            hStr.clear();
        } else {
            // hStr becomes the remaining heuristics in the chain
            hStr = hStr.substr(hyphenLoc + 1);
        }
    }
    // Return the newly populated heuristic chain for easy reference
    return true;
}

void
HeuristicChain::setESTList(ESTList* estList) {
    this->estList = estList;
}

int
HeuristicChain::run() {
    int retVal = 0;
    for (size_t i = 0; (i < chain.size()); i++) {
	if ((retVal = chain[i]->run()) != 0) {
	    // Immediately stop when a heuristic fails to run successfully
	    break;
	}
    }
    // Return zero if all heuristics ran successfully.
    return retVal;
}

#endif
