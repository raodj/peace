#ifndef INPUT_FILE_CPP
#define INPUT_FILE_CPP

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

#include "InputFile.h"
#include <cstdlib>

InputFile::InputFile(const std::string& srcFileName) :
    fileName(srcFileName) {
    randomizeNbases = false;
    maskBases       = false;
    normalizeNTs    = true;
}

InputFile::~InputFile() {
    // Nothing to be done for now.
}

EST*
InputFile::getNextEntry(const int estID) {
    // First check to ensure that we currently don't have this estID
    // already in our list.
    if (entryOffset.find(estID) != entryOffset.end()) {
        // The entry already exists. Cannot add duplicate entry
        return NULL;
    }
    // Save the location from where we are reading data for this entry
    const size_t position = getCurrPos();
    // Read information about the fragment using helper method.
    std::string info, sequence;
    QualityVector quality;
    if (!readEntry(info, sequence, quality)) {
        // Error occured when reading the sequence.
        return NULL;
    }
    // Normalize the sequence by applying makign and ranodmization
    // options via a helper method.
    if (normalizeNTs) {
        sequence = EST::normalizeBases(sequence, maskBases,
                                       randomizeNbases, (int) position);
    }
    // OK, create and return the EST after updating our internal tables
    entryOffset[estID] = position;
    // Normalize the nucleotide sqeuence 
    return new EST(estID, info, sequence, quality);
}

bool
InputFile::repopulate(EST& est) {
    // First find the position where the entry occurs in the input file.
    std::map<int, size_t>::const_iterator entry = entryOffset.find(est.getID());
    if (entry == entryOffset.end()) {
        // The entry does not exists. Cannot repopulate
        return false;
    }
    // Have the derived class implementation seek to appropriate position
    if (!setCurrPos(entry->second)) {
        // Seeking to the previously saved position failed!
        return false;
    }
    // Now, have the helper method reload the core information
    std::string info, sequence;
    QualityVector quality;
    if (!readEntry(info, sequence, quality)) {
        // Error occured when reading the sequence.
        return NULL;
    }
    // Normalize the sequence by applying makign and ranodmization
    // options via a helper method.
    sequence = EST::normalizeBases(sequence, maskBases,
                                   randomizeNbases, (int) entry->second);
    // Update the information in the EST
    est.repopulate(info, sequence, quality);
    // All the processing was successfully completed
    return true;
}

void
InputFile::setOptions(bool maskBases, bool randomizeNbases, bool normalizeNTs) {
    this->maskBases       = maskBases;
    this->randomizeNbases = randomizeNbases;
    this->normalizeNTs    = normalizeNTs;
}

#endif
