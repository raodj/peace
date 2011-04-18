#ifndef PARAMETER_SET_MANAGER_CPP
#define PARAMETER_SET_MANAGER_CPP

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

#include "ParameterSetManager.h"
#include "Utilities.h"
#include "ESTList.h"

#include <cstring>
#include <climits>

ParameterSetManager::ParameterSetManager() {
    // Initialize pointers to default initial value
    estList     = NULL;
    lookupTable = NULL;
}

void
ParameterSetManager::setupParameters(int t1, int u1, int ws1, int t2,
                                     int u2, int ws2, int t3, int u3,
                                     int ws3) {
    // Clear out any earlier values
    resetParameterSets();
    // Add the default parameter sets to the list of sets.
    // Line 1 is Two-Pass D2, Line 2 is heuristic parameters
    addParameterSet(new ParameterSet(-1, 150, 50, 1, 40, 40,
                                     t1, u1, ws1));
    addParameterSet(new ParameterSet(150, 400, 75, 25, 60, 100,
                                     t2, u2, ws2));
    addParameterSet(new ParameterSet(400, INT_MAX, 100, 50, 80,
                                     130, t3, u3, ws3));
}

void
ParameterSetManager::setupParameters(const bool) {
    // Clear out any earlier values
    resetParameterSets();
    // Add the default parameter set to the list of sets.
    addParameterSet(new ParameterSet(-1, INT_MAX, 100, 50, 40,
                                     130, 65, 6, 8));
}

void
ParameterSetManager::addParameterSet(ParameterSet* p) {
    parameterSets.push_back(p);
}

ParameterSetManager::~ParameterSetManager() {
    // Clear out look up table entries.
    resetLookupTable();
    resetParameterSets();
}

void
ParameterSetManager::resetParameterSets() {
    // Now free up parameter sets
    for(size_t i = 0; (i < parameterSets.size()); i++) {
        delete parameterSets[i];
    }
    parameterSets.clear();
}

void
ParameterSetManager::resetLookupTable() {
    // Clear out our look-up tables as needed
    if (lookupTable != NULL) {
        for(size_t i = 0; (i <= parameterSets.size()); i++) {
            delete [] lookupTable[i];
        }
        delete [] lookupTable;
        lookupTable = NULL;
    }
}

void
ParameterSetManager::initialize() {
    // Reset any earlier values lookup table entries
    resetLookupTable();
    // If we don't have any parameter sets added, then we setup
    // default values here.
    if (parameterSets.empty()) {
        setupParameters(false);
    }
    // Note that this method may be called multiple times!
    // First populate the lengthToParamSetMap vector.
    ASSERT( estList != NULL );
    const int estCount = estList->size();
    // Reserve sufficient space to minimize overheads
    lengthToParamSetMap.reserve(estCount);
    // Now compute the mappings using a helper method.
    for(int index = 0; (index < estCount); index++) {
        const EST* est   = estList->get(index);
        const int estLen = est->getSequenceLength();
        lengthToParamSetMap.push_back(getParameterSetIndex(estLen));
    }
    // Now create and populate the look up tables based on the entries
    // in the parameterSets.  First create the dynamic 2-D array.
    const int setCount = parameterSets.size();
    lookupTable = new int*[setCount + 1];
    for(int row = 0; (row <= setCount); row++) {
        lookupTable[row] = new int[setCount + 1];
    }
    // Populate the entries.  Note that the last row and column
    // correspond to invalid entries and they are treated as a
    // special/edge case.
    for(int index1 = 0; (index1 < setCount + 1); index1++) {
        for(int index2 = 0; (index2 < setCount + 1); index2++) {
            if (std::max(index1, index2) == setCount) {
                // This is an edge case that corresponds to invalid
                // entries maintained as the last row and column.
                lookupTable[index1][index2] = -1;
            } else {
                lookupTable[index1][index2] =
                    getPreferredSetIndex(index1, index2);
            }
        }
    }
}

int
ParameterSetManager::getPreferredSetIndex(const int index1,
                                          const int index2) const {
    int prefIndex = index1;
    if (parameterSets[prefIndex]->minLength >
        parameterSets[index2]->minLength) {
        prefIndex = index2;
    }
    // Check to ensure we are not comparing the shortest and longest
    // entries as it does not make much sense.
    const int setCount = parameterSets.size();
    if ((prefIndex == 0) && (setCount > 2) &&
        (std::max(index1, index2) == setCount - 1)) {
        // Don't compare these two entries.
        return -1;
    }
    // Return the preferred index of the two.
    return prefIndex;
}

int
ParameterSetManager::getParameterSetIndex(const int estLength) const {
    const int setCount = (int) parameterSets.size();
    for (int setNum = 0; (setNum < setCount); setNum++) {
        // We assume that the parameter sets are already organized in
        // increasing order of fragment lengths.
        if (estLength <= parameterSets[setNum]->maxLength) {
            return setNum;
        }
    }
    
    return setCount;
}

void
ParameterSetManager::setESTList(const ESTList* list) {
    estList = list;
}

void
ParameterSetManager::entriesAdded (int startIndex, int endIndex) {
    // We currently have a rather simple implementation here...
    ASSERT ( startIndex == (int) lengthToParamSetMap.size() );
    for(int i = startIndex; (i < endIndex); i++) {
        const EST* est   = estList->get(i);
        ASSERT( est != NULL );
        const int estLen = est->getSequenceLength();
        // Add new entry to our look-up table.
        lengthToParamSetMap.push_back(getParameterSetIndex(estLen));
    }
}

void
ParameterSetManager::entriesRemoved(int startIndex, int endIndex) {
    // We currently have a rather simple implementation here...    
    ASSERT( endIndex == (int) lengthToParamSetMap.size() );
    for(int i = endIndex; (i > startIndex); i--) {
        lengthToParamSetMap.pop_back();
    }
}

#endif
