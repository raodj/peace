#ifndef EST_LIST_CPP
#define EST_LIST_CPP

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

#include "ESTList.h"

ESTList::ESTList() {
    // All other members initialize themselves.
    maxESTLen = 0;
}

ESTList::~ESTList() {
    // Clear out all the information we currently have
    ESTList::reset();
}

EST*
ESTList::add(const int id, const std::string& info,
             const std::string& sequence) {
    QualityVector quality;
    // Use overloaded method to do the actual job
    return add(id, info, sequence, quality);
}

EST*
ESTList::add(const int id, const std::string& info,
             const std::string& sequence, const QualityVector& quality) {
    if (id != (int) estVector.size()) {
        // This id is not acceptable. Sorry.
        return NULL;
    }
    // Update new max EST length (if it has been counted)
    if (maxESTLen != 0) {
        // Update max length
        maxESTLen = std::max(maxESTLen, (int) sequence.size());
    }
    // Instantiate new EST
    EST *newEST = new EST(id, info, sequence, quality);
    // Add est to end of the est list.
    estVector.push_back(newEST);
    // return the newly created EST back tot he caller
    return newEST;
}

void
ESTList::reset() {
    const int ESTCount = (int) estVector.size();
    for(int id = 0; (id < ESTCount); id++) {
        delete estVector[id];
    }
    // Clear out all entries as they have all been deleted.
    estVector.clear();
    // Reset other variables & counters
    maxESTLen = 0;
}

void
ESTList::deleteLastESTs(const int count) {
    for(int i = 0; ((estVector.size() > 0) && (i < count)); i++) {
        // Remove the last entry
        delete estVector.back();
        estVector.pop_back();
    }
}

void
ESTList::dumpESTList(std::ostream& os) const {
    const int EstCount = (int) estVector.size();
    for(int id = 0; (id < EstCount); id++) {
        // Dump the EST information out.
        estVector[id]->dumpEST(os);
    }
}

void
ESTList::dumpESTList(std::ostream& os, const bool processed) const {
    const int EstCount = (int) estVector.size();
    for(int id = 0; (id < EstCount); id++) {
        // Dump the EST information out based on processed status.
        if (estVector[id]->hasBeenProcessed() == processed) {
            // Matched condition.
            estVector[id]->dumpEST(os);
        }
    }
}

int
ESTList::getProcessedESTCount() const {
    int count = 0;
    const int ESTCount = (int) estVector.size();
    for(int id = 0; (id < ESTCount); id++) {
        if (estVector[id]->hasBeenProcessed()) {
            count++;
        }
    }
    return count;
}

size_t
ESTList::getMaxESTLen() const {
    if (maxESTLen == 0) {
        // First time we have to compute the max length as it is not
        // yet known.
        const int ESTCount = (int) estVector.size();
        for(int id = 0; (id < ESTCount); id++) {
            maxESTLen = std::max(maxESTLen, estVector[id]->getSequenceLength());
        }
    }
    // We have the maximum length EST
    return maxESTLen;
}

const EST*
ESTList::repopulate(int index) const {
    return estVector[index];
}

bool
ESTList::add(InputFile* UNREFERENCED_PARAMETER(inputFile),
             const long UNREFERENCED_PARAMETER(startIndex),
             const long UNREFERENCED_PARAMETER(endIndex)) {
    return false;
}

#endif
