#ifndef ON_DEMAND_EST_LIST_CPP
#define ON_DEMAND_EST_LIST_CPP

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

#include "OnDemandESTList.h"
#include "InputFile.h"
#include "InputFileFactory.h"

OnDemandESTList::OnDemandESTList() {
    // All other members initialize themselves.
}

OnDemandESTList::~OnDemandESTList() {
    // Clear out all the information we currently have
    reset();
}

bool
OnDemandESTList::add(InputFile* inputFile,
                     const long startIndex, const long endIndex) {
    // Read all the entries from the input file and add it to our EST
    // list.
    const int ESTStartID = estVector.size();
    for(int index = 0; ((inputFile->good()) && inputFile->hasNextEntry());
        index++) {
        // Read the next entry from the file
        EST *est = inputFile->getNextEntry(index + ESTStartID);
        if (est == NULL) {
            // Huh. error occured when loading EST. Exit ASAP.
            std::cerr << "Error reading cDNA fragment from data file: "
                      << inputFile->getFileName() << std::endl;
            return false;
        }
        // Unpopulate entries if needed
        if ((est->getID() < startIndex) || (est->getID() > endIndex)) {
            est->unpopulate();
        }
        // Add it to our vector.
        estVector.push_back(est);
    }
    // Return true if the file is still good at the end of processing
    // all the entries in the input file.
    if (inputFile->good()) {
        inputFileList.push_back(inputFile);
        return true;
    }
    // Something went wrong when reading the inputs
    return false;
}

EST*
OnDemandESTList::repopulate(int estID) const {
    if ((estID < 0) || (estID >= (int) estVector.size())) {
        // Invalid EST id
        return NULL;
    }
    // Extract the entry.
    EST *est = estVector[estID];
    // Search through the files trying to repopulate the EST.
    for(size_t id = 0; (id < inputFileList.size()); id++) {
        if (inputFileList[id]->repopulate(*est)) {
            // Successfully repopulated
            return est;
        }
    }
    // Could not repopulate the entry
    return NULL;
}

void
OnDemandESTList::reset() {
    // Close any input files we have and clear those entries
    for(size_t index = 0; (index < inputFileList.size()); index++) {
        delete inputFileList[index];
    }
    inputFileList.clear();
    // Let base class to its task
    ESTList::reset();
}

std::string
OnDemandESTList::getESTInfo(const int estID) const {
    if ((estID < 0) || (estID >= (int) estVector.size())) {
        // Invalid EST id
        return "";
    }
    // Extract the entry.
    const EST *est = estVector[estID];
    // Check if we already have the data
    if (est->isPopulated()) {
        return est->getInfo();
    }
    // Create a mutable copy of the requested EST
    EST copy(*est);
    // Search through the files trying to repopulate the EST.
    for(size_t id = 0; (id < inputFileList.size()); id++) {
        if (inputFileList[id]->repopulate(copy)) {
            // Successfully repopulated
            return copy.getInfo();
        }
    }
    // Could not repopulate the entry
    return "";
}

#endif
