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
#include "FilterFactory.h"
#include "Utilities.h"
#include "MPIHelper.h"
#include "MPIStats.h"
#include "EST.h"

// Get rid of magic numbers
#define NO_ERROR 0

// The static pointer to the singleton filter chain instance.
FilterChain* FilterChain::ptrInstance = NULL;

void
FilterChain::showArguments(std::ostream& os) {
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->showArguments(os);
    }
}

bool
FilterChain::parseArguments(int& argc, char **argv) {
    bool flag = true;
    for (size_t i = 0; (i < chain.size()); i++) {
        if (!chain[i]->parseArguments(argc, argv)) {
            flag = false;
        }
    }
    return flag;
}

bool
FilterChain::addFilter(Filter* h) {
    ASSERT( h != NULL );
    chain.push_back(h);
    // Everything went well
    return true;
}

int
FilterChain::initialize() {
    int result;
    for (size_t i = 0; (i < chain.size()); i++) {
        if ((result = chain[i]->initialize()) != NO_ERROR) {
            // Error occured during initialization. Bail out.
            std::cerr << "Error initializing filter "
                      << chain[i]->getName() << std::endl;
            return result;
        }
    }
    return 0; // Everything went well
}

void
FilterChain::finalize() {
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->finalize();
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

FilterChain::~FilterChain() {
    for(size_t i = 0; (i < chain.size()); i++) {
        delete chain[i];
    }
    chain.clear();
}

FilterChain::FilterChain() {
    // Nothing to be done for now
}

void
FilterChain::printStats(std::ostream& os, const int rank) const {
    os << "Filters on process with Rank " << rank << "\n"
       << "-------------------------------------------\n";
    for (size_t i = 0; (i < chain.size()); i++) {
        chain[i]->printStats(os);
    }
}


FilterChain*
FilterChain::setupChain(const char* filterStr, ClusterMaker *clusterMaker) {
    // First validate and create a blank filter chain.
    ASSERT ( ptrInstance == NULL );
    ptrInstance = new FilterChain();
    // Check if there is anything to be created.
    if (filterStr == NULL) {
        return NULL;
    }
    // Process the non-empty filter string.
    std::string fStr(filterStr);
    FilterChain* filterChain = ptrInstance;
    ASSERT ( filterChain != NULL );
    // Process one word at a time. Words are separated by a "hyphen"
    // character.
    while (!fStr.empty()) {
        // Locate the next hyphen character and get name of filter
        const std::string::size_type hyphenLoc = fStr.find('-');
        const std::string name = fStr.substr(0, hyphenLoc);
        // Create filter and add it to the chain
        Filter *filter =
            FilterFactory::create(name.c_str(), clusterMaker);
        if (filter == NULL) {
            // Break out and return null, which will cause main to show usage
            return NULL;
        }
        // This assert caused errors if invalid filters were entered
        //ASSERT( filter != NULL );
        filterChain->addFilter(filter);
        // Remove already created filter to process next one in chain.
        if (hyphenLoc == std::string::npos) {
            // No hyphen, so this was the last filter in the chain
            fStr.clear();
        } else {
            // fStr becomes the remaining filters in the chain
            fStr = fStr.substr(hyphenLoc + 1);
        }
    }
    // Return the newly populated filter chain for easy reference
    return filterChain;
}

void
FilterChain::getOwnedESTidx(int& startIndex, int& endIndex) {
    // Exact implementation as in
    // MSTClusterMaker::getOwnedESTidx. This method has been
    // copy-pasted so that filters can operate on their own different
    // sub-set of ESTs if they choose. Maybe the method can be
    // combined together.
    const int ESTsPerProcess = EST::getESTList().size() / MPI_GET_SIZE();
    const int ExtraESTs      = EST::getESTList().size() % MPI_GET_SIZE();
    const int MyRank         = MPI_GET_RANK();
    
    // First figure out the starting and ending EST this processs is
    // responsible for further use.
    startIndex = MyRank * ESTsPerProcess;
    // Check for extra proceses as needed.
    if (MyRank <= ExtraESTs) {
        // The previous processes have one extra ESTs as number of
        // ESTs are not evenly divisible by number of processes.  So
        // account for this.
        startIndex = ((ESTsPerProcess + 1) * MyRank);
    } else {
        startIndex += ExtraESTs;
    }
    
    // Compute the last est index this process owns.
    endIndex = startIndex + ESTsPerProcess;
    if (MyRank < ExtraESTs) {
        // This process owns one extra EST as ESTs are not evenly
        // divisible by the number of processes.
        endIndex++;
    }
}

int
FilterChain::applyFilters(ClusterMaker *clusterMaker) {
    ASSERT ( ptrInstance != NULL );
    // Save original number of ESTs for sanity check below.
    const int prevEstCount = EST::getESTCount();
    // First initialize all the filters.
    int resultCode;
    if ((resultCode = ptrInstance->initialize()) != 0) {
        // Error occured bail out.
        return resultCode;
    }
    // Apply filters on the subset of ESTs that we are responsible for
    // on this process.
    int startIndex, endIndex;
    getOwnedESTidx(startIndex, endIndex);
    for(int estIdx = startIndex; (estIdx < endIndex); estIdx++) {
        // User helper method to apply the filters.
        if (!EST::getEST(estIdx)->hasBeenProcessed()) {
            ptrInstance->applyFilters(estIdx);
        }
    }
    // Now participate in global iterative broadcast chain. Possibly
    // this loop can be replaced by a MPI all-to-all broadcast that
    // expense of increased memory footprint.
    if (MPI_GET_SIZE() > 1) {
        allToAllBroadcast(clusterMaker);
    }    
    // Finalize all the filters
    ptrInstance->finalize();
    // Ensure numebr of ESTs are back to normal
    if (prevEstCount != EST::getESTCount()) {
        std::cerr << "The filters did not remove all the dummy ESTs they "
                  << "added. This is a programming error. Aborting!\n";
        resultCode = 2;
    }
    
    // All operations proceeded successfully
    return resultCode;
}

void
FilterChain::allToAllBroadcast(ClusterMaker *clusterMaker) {
    const int ProcCount = MPI_GET_SIZE();
    ASSERT (ProcCount > 1);
    // First build our local filter information for broadcast when it
    // is our turn.
    std::vector<int> localSuperList;
    for(size_t i = 0; (i < ptrInstance->chain.size()); i++) {
        ptrInstance->chain[i]->addFilterData(localSuperList);
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
