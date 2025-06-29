#ifndef PDB_LIST_FILE_CPP
#define PDB_LIST_FILE_CPP

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

#include <sstream>
#include "PDBListFile.h"

PDBListFile::PDBListFile(const std::string& fileName) : InputFile(fileName) {
    // Open the FASTA file to be read and processed
    pdbListFile.open(fileName.c_str());
}

PDBListFile::~PDBListFile() {
    pdbListFile.close();
}

bool
PDBListFile::hasNextEntry() {
    return !pdbListFile.eof() && (pdbListFile.peek() != EOF);
}

size_t
PDBListFile::getCurrPos() {
    return pdbListFile.tellg();
}

bool
PDBListFile::setCurrPos(const size_t position) {
    pdbListFile.seekg(position, std::ios::beg);
    if (!pdbListFile.fail()) {
        // The seek operation was successful. If the stream was not in
        // a good state earlier (because we hit EOF), we reset the EOF
        // flag here to permit smooth reading after this seek.
        pdbListFile.clear();
        return true;
    }
    // Error during seek
    return false;
}

bool
PDBListFile::readEntry(std::string& info, std::string& sequence,
                       QualityVector& UNUSED(quality)) {
    // Some basic validation on the fastFile first.    
    if (!pdbListFile.good()) {
        // Can't read information from the fasta file.
        return false;
    }

    // Skip over any empty and comment lines to get to a valid
    // non-comment line.
    std::string line;
    while (std::getline(pdbListFile, line)) {
        trim(line);
        if (!line.empty() && (line[0] != '#')) {
            // This is a valid PDB file entry. Use it.
            info = line;
            // Returns true if PDB data was successfully loaded into
            // sequence.
            return loadPDB(line, sequence);
        }
    }
    // Return false to indicate a valid entry was not found.
    return false;
}


bool
PDBListFile::loadPDB(const std::string& path, std::string& data) const {
    // Open the specified PDB file and ensure it went well.
    std::ifstream pdb(path);
    if (!pdb.good()) {
        return false;  // Error opening pdb file.
    }
    // Load the whole file into data.
    std::stringstream buffer;
    buffer << pdb.rdbuf();
    data = buffer.str();
    return !data.empty();  // True if we got some data.
}

#endif
