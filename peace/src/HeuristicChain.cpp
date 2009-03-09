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
// Authors: James C. Moler         molerjc@muohio.edu
//
//---------------------------------------------------------------------------

#include "HeuristicChain.h"

// The static list of heuristics in the heuristic chain.
std::vector<Heuristic*> HeuristicChain::chain;

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

bool
HeuristicChain::addHeuristic(Heuristic* h) {
    ASSERT( h != NULL );
    chain.push_back(h);
    // Everything went well
    return true;
}

bool
HeuristicChain::shouldAnalyze() {
    // Return true as a default
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

#endif
