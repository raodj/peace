#ifndef FILTER_CPP
#define FILTER_CPP

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

#include "EST.h"
#include "Filter.h"
#include "Utilities.h"
#include "ClusterMaker.h"
#include "RuntimeContext.h"

Filter::Filter(const std::string& name, RuntimeContext* context)
    : filterName(name), runtimeContext(context) {
    runCount    = 0;
    filterCount = 0;
    // Nothing else to be done for now.
}

Filter::~Filter() {
    // Empty constructor begets an empty destructor
}

bool
Filter::applyFilter(const EST* est) {
    int clusterIndex = -1;
    bool hasBeenFiltered = false;
    // Check if the derived filter class tells us if this EST must be
    // filtered out.
    if ((hasBeenFiltered = runFilter(est, clusterIndex))) {
        // This EST was filtered out.
        filterCount++;
        // Track the filtered EST information for distribution to
        // other processes later on.
        filteredESTList[clusterIndex].push_back(est->getID());
        // Add this EST to the clsuter (in cluster maker)
        ClusterMaker *clusterMaker = runtimeContext->getClusterMaker();
        if (clusterMaker != NULL) {
            clusterMaker->addEST(clusterIndex, est->getID());
        }
    }
    // Track number of times this filter was run.
    runCount++;
    // Return cluster index or -1 if this EST passed this filter test.
    return hasBeenFiltered;
}

void
Filter::addFilterData(std::vector<int>& superList) const {
    for(FilteredESTList::const_iterator curr = filteredESTList.begin();
        (curr != filteredESTList.end()); curr++) {
        // Add cluster index to super list first.
        superList.push_back(curr->first);
        // Now add information about the ests 
        const std::vector<int>& estList = curr->second;
        // Add number of ests filtered
        superList.push_back(estList.size());
        // Add the ESTs themselves.
        superList.reserve(estList.size());
        superList.insert(superList.end(), estList.begin(), estList.end());
    }
}

void
Filter::processFilterData(const std::vector<int>& superList,
                          ClusterMaker *clusterMaker) {
    // Process sets of entries until we find a end-of-list marker with
    // the entry "-1, 0".
    int currIndex = 0;
    while (superList[currIndex] != -1) {
        // First entry is the cluster to which ESTs are to be added.
        const int clusterID = superList[currIndex];
        // The second entry is number of ESTs in the list.
        const int estCount  = superList[currIndex + 1];
        // Add all ESTs to the dummy cluster (if cluster maker is
        // specified).
        if (clusterMaker != NULL) {
            for(int i = 0; (i < estCount); i++) {
                clusterMaker->addEST(clusterID, superList[currIndex + i + 2]);
            }
        }
        // We have processed one batch of filter data. Onto the next.
        currIndex += (estCount + 2);
    }
}

void
Filter::printStats(std::ostream& os) const {
    os << "Statistics from " << getName()   << " filter:\n"
       << "\tNumber of calls         : " << runCount     << "\n"
       << "\tNumber of ESTs filtered : " << filterCount  << std::endl;
}

int
Filter::addDummyCluster(const std::string& clusterInfo) {
    ASSERT( runtimeContext != NULL );
    ClusterMaker *clusterMaker = runtimeContext->getClusterMaker();
    if ( clusterMaker == NULL ) {
        // No cluster maker available to add a dummy cluster
        return -1;
    }
    // Add dummy cluster to the cluster maker
    return clusterMaker->addDummyCluster(clusterInfo);
}

#endif
