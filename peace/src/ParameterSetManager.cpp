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

// The static pointer to the singleton parameter set manager instance.
ParameterSetManager* ParameterSetManager::ptrInstance = NULL;

void
ParameterSetManager::setupParameters(int t1, int u1, int ws1, int t2,
                                     int u2, int ws2, int t3, int u3,
                                     int ws3) {
    // First validate and create a blank parameter set manager.
    ASSERT ( ptrInstance == NULL );
    ptrInstance = new ParameterSetManager;

    // Add the default parameter sets to the list of sets.
    // Line 1 is Two-Pass D2, Line 2 is heuristic parameters
    ptrInstance->addParameterSet(new ParameterSet(-1, 150, 50, 1, 45, 45,
                                                  t1, u1, ws1));
    ptrInstance->addParameterSet(new ParameterSet(150, 400, 75, 25, 75, 100,
                                                  t2, u2, ws2));
    ptrInstance->addParameterSet(new ParameterSet(400, -1, 100, 50, 105, 130,
                                                  t3, u3, ws3));
}

void
ParameterSetManager::addParameterSet(ParameterSet* p) {
    parameterSets.push_back(p);
}

ParameterSetManager::ParameterSetManager() { 
    // Nothing more to do
}

ParameterSetManager::~ParameterSetManager() {
    for(size_t i = 0; (i < parameterSets.size()); i++) {
        delete parameterSets[i];
    }
    parameterSets.clear();
}

ParameterSet*
ParameterSetManager:: getParameterSet(const int seq1Len, const int seq2Len) {
    int minLength = std::min(seq1Len, seq2Len);
    int setNum = 0;
    for (setNum = 0; setNum < parameterSets.size()-1; setNum++) {
        if (minLength <= parameterSets[setNum]->maxLength) {
            break;
        }
    }
    if (setNum == 0 && parameterSets.size() > 2 &&
        std::max(seq1Len, seq2Len) >
        parameterSets[parameterSets.size()-1]->minLength) {
        // Sequences belong in the least and greatest parameter sets
        // Therefore, no comparison should be made
        return NULL;
    } else {
        return parameterSets[setNum];
    }
}

int
ParameterSetManager::getMaxFrameSize() {
    return parameterSets[parameterSets.size()-1]->frameSize;
}

#endif
