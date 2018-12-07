#ifndef FILTER_CHAIN_CPP
#define FILTER_CHAIN_CPP

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

#include "FilterChain.h"
#include "Utilities.h"
#include "MPIHelper.h"
#include "MPIStats.h"
#include "ESTList.h"
#include "SubSystem.h"
#include "RuntimeContext.h"
#include <cstdlib>

FilterChain::FilterChain() : Component("FilterChain") {
    origESTCount = 0;
}


FilterChain::~FilterChain() {
    for(size_t i = 0; (i < chain.size()); i++) {
        delete chain[i];
    }
    chain.clear();
}

bool
FilterChain::addFilter(Filter* h) {
    ASSERT( h != NULL );
    chain.push_back(h);
    // Everything went well
    return true;
}

bool
FilterChain::initialize() {
    // Avoid redundant initialization.
    if (isInitialized()) {
        // Already initialized. Nothing further to be done.
        return true;
    }
    // Initialize the base class first
    if (!Component::initialize()) {
        return false;
    }
    // Save original number of ESTs for sanity check in finalize() method
    ASSERT ( subSystem != NULL );
    ASSERT ( subSystem->getContext() != NULL );
    ASSERT ( subSystem->getContext()->getESTList() != NULL );        
    origESTCount = subSystem->getContext()->getESTList()->size();

    // Initialize all the filters in the chain.
    int result;
    for (size_t i = 0; (i < chain.size()); i++) {
        if ((result = chain[i]->initialize()) != 0) {
            // Error occured during initialization. Bail out.
            std::cerr << "Error initializing filter "
                      << chain[i]->getName() << std::endl;
            return false;
        }
    }
    return true; // Everything went well
}

int
FilterChain::run() {
    ASSERT ( subSystem != NULL );
    ASSERT ( subSystem->getContext() != NULL );
    ASSERT ( subSystem->getContext()->getESTList() != NULL );    
    // Get convenience reference to global list of ESTs
    ESTList& estList = *(subSystem->getContext()->getESTList());
    // Apply filters on the subset of ESTs that we are responsible for
    // on this process.
    int startIndex, endIndex;
    getOwnedESTidx(estList.size(), startIndex, endIndex);
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        // Use helper method to apply the filters.
        EST *est = estList.get(estIdx, true);
        if (!est->hasBeenProcessed()) {
            // Filters track ESTs they have filtered out.
            applyFilters(est);
	    estList.unpopulate(estIdx);
        }
    }
    // Now participate in global iterative broadcast chain. Possibly
    // this loop can be replaced by a MPI all-to-all broadcast that
    // expense of increased memory footprint.
    if (MPI_GET_SIZE() > 1) {
        ClusterMaker *clusterMaker = subSystem->getContext()->getClusterMaker();
        allToAllBroadcast(clusterMaker);
    }
    // All operations proceeded successfully
    return 0;
}

void
FilterChain::finalize() {
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->finalize();
    }
    // Ensure numebr of ESTs are back to normal
    ASSERT ( subSystem != NULL );
    ASSERT ( subSystem->getContext() != NULL );
    ASSERT ( subSystem->getContext()->getESTList() != NULL );    
    // Get convenience reference to shared list of ESTs
    const ESTList* estList = subSystem->getContext()->getESTList();
    if (origESTCount != estList->size()) {
        std::cerr << "The filters did not remove all the dummy ESTs they "
                  << "added. This is a programming error. Aborting!\n";
        exit(3);
    }
}

Filter*
FilterChain::getFilter(const std::string& name) const {
    for(size_t i = 0; (i < chain.size()); i++) {
        if (chain[i]->getName() == name) {
            // Found a match!
            return chain[i];
        }
    }
    // A matching filter was not found.
    return NULL;
}

void
FilterChain::printStats(std::ostream& os, const int rank) const {
    os << "Filters on process with Rank " << rank << "\n"
       << "-------------------------------------------\n";
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->printStats(os);
    }
}

void
FilterChain::allToAllBroadcast(ClusterMaker *clusterMaker) {
    const int ProcCount = MPI_GET_SIZE();
    ASSERT (ProcCount > 1);
    // First build our local filter information for broadcast when it
    // is our turn.
    std::vector<int> localSuperList;
    for(size_t i = 0; (i < chain.size()); i++) {
        chain[i]->addFilterData(localSuperList);
    }
    // Add list termination entries.
    localSuperList.push_back(-1);
    localSuperList.push_back(0);
    
    // Receive broadcast from all other processes. But
    // when it is our turn we broadcast our local filter information.
    for(int rank = 0; (rank < ProcCount); rank++) {
        if (rank == MPI_GET_RANK()) {
            // It is our turn to broadcast the information.
            int size = localSuperList.size();
            ASSERT( size > 0 );
            MPI_BCAST(&size, 1, MPI_INT, rank);
            MPI_BCAST(&localSuperList[0], size, MPI_INT, rank); 
        } else {
            // Receive from another process and merge filter data.
            int size = 0;
            MPI_BCAST(&size, 1, MPI_INT, rank);
            ASSERT ( size > 0 );
            std::vector<int> remoteSuperList;
            remoteSuperList.reserve(size);
            MPI_BCAST(&remoteSuperList[0], size, MPI_INT, rank);
            // Process the remote list
            Filter::processFilterData(remoteSuperList, clusterMaker);
        }
    }
}

#endif
