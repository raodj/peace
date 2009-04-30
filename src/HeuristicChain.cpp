#ifndef HEURISTIC_CHAIN_CPP
#define HEURISTIC_CHAIN_CPP

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
//          James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "HeuristicChain.h"
#include "HeuristicFactory.h"
#include "Utilities.h"

// Get rid of magic numbers
#define NO_ERROR 0

// The static pointer to the singleton heuristic chain instance.
HeuristicChain* HeuristicChain::ptrInstance = NULL;

HeuristicChain*
HeuristicChain::getHeuristicChain() {
    if (ptrInstance == NULL) {
        // Initialize the heuristic chain
        ptrInstance = new HeuristicChain;
    }
    ASSERT ( ptrInstance != NULL );
    return ptrInstance;
}

void
HeuristicChain::showArguments(std::ostream& os) {
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->showArguments(os);
    }   
}

bool
HeuristicChain::parseArguments(int& argc, char **argv) {
    bool flag = true;
    for (size_t i = 0; (i < chain.size()); i++) {
        if (!chain[i]->parseArguments(argc, argv)) {
            flag = false;
        }
    }
    return flag;
}

bool
HeuristicChain::addHeuristic(Heuristic* h) {
    ASSERT( h != NULL );
    chain.push_back(h);
    // Everything went well
    return true;
}

int
HeuristicChain::initialize() {
    int result;
    for (size_t i = 0; (i < chain.size()); i++) {
        if ((result = chain[i]->initialize()) != NO_ERROR) {
            // Error occured during initialization. Bail out.
            std::cerr << "Error initializing heuristic "
                      << chain[i]->getName() << std::endl;
            return result;
        }
    }
    return 0; // Everything went well
}

int
HeuristicChain::setReferenceEST(const int estIdx) {
    int result;
    for (size_t i = 0; (i < chain.size()); i++) {
        if ((result = chain[i]->setReferenceEST(estIdx)) != NO_ERROR) {
            // Error occured during initialization. Bail out.
            return result;
        }
    }
    return 0; // Everything went well
}

bool
HeuristicChain::shouldAnalyze(const int otherEST) {
    for (size_t i = 0; (i < chain.size()); i++) {
        if (!chain[i]->shouldAnalyze(otherEST)) {
            // Immediately stop when a heuristic says no
            return false;
        }
    }
    // If all heuristics say yes, return true
    return true;
}

HeuristicChain::~HeuristicChain() {
    for(size_t i = 0; (i < chain.size()); i++) {
        delete chain[i];
    }
    chain.clear();
}

HeuristicChain::HeuristicChain() {
    // Nothing to be done for now
}

HeuristicChain*
HeuristicChain::setupChain(const char* heuristicStr, const int refESTidx,
                           const std::string& outputFile) {
    if (heuristicStr == NULL) {
        return NULL;
    }
    // Process the non-empty heuristic string.
    std::string hStr(heuristicStr);
    HeuristicChain* heuristicChain = HeuristicChain::getHeuristicChain();
    ASSERT ( heuristicChain != NULL );
    // Process one word at a time. Words are separated by a "hyphen"
    // character.
    while (!hStr.empty()) {
        // Locate the next hyphen character and get name of heuristic
        const std::string::size_type hyphenLoc = hStr.find('-');
        const std::string name = hStr.substr(0, hyphenLoc);
        // Create heuristic and add it to the chain
        Heuristic *heuristic =
            HeuristicFactory::create(name.c_str(), refESTidx, outputFile);
        ASSERT( heuristic != NULL );
        heuristicChain->addHeuristic(heuristic);
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
    return heuristicChain;
}

#endif
