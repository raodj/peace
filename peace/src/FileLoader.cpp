#ifndef FILE_LOADER_CPP
#define FILE_LOADER_CPP

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
// as a result of using, modifying or distributing this software or
// its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "FileLoader.h"
#include "SFFReader.h"
#include "MPIHelper.h"
#include "EST.h"

// The sentinel file name that indicates not to load a data file
const std::string FileLoader::IgnoreFileName("<none>");

bool
FileLoader::loadFASTAFile(const std::string& toolName, const char *fileName,
                          const bool maskBases, const bool randomizeNbases) {
    if (FileLoader::IgnoreFileName == fileName) {
        // The user does not want to use a file name. This possibly
        // happens when using peace libary from custom applications.
        return true;
    }
    FILE *fastaFile = NULL;
#ifndef _WINDOWS
    fastaFile = fopen(fileName, "rt");
#else
    fopen_s(&fastaFile, fileName, "rt");
#endif
    if ((fastaFile == NULL) || (ferror(fastaFile))) {
        std::cerr << toolName << "(Rank: ";
        std::cerr << MPI_GET_RANK()
                  << "): Error opening FASTA file "
                  << fileName << " for reading." << std::endl;
        return false;
    }
    // Track line number in file being processed.
    int lineNum = 1;
    // Repeatedly read EST's from the file.
    while (!feof(fastaFile)) {
        EST *est = EST::create(fastaFile, lineNum, maskBases,
                               randomizeNbases);
        if ((est == NULL) && (!feof(fastaFile) || ferror(fastaFile))) {
            // An error occured when reading EST.
            fclose(fastaFile);
            std::cerr << toolName << ": Error loading EST from "
                      << "FASTA file " << fileName << " at line: "
                      << lineNum << std::endl;
            return false;
        }
    }
    // All the EST's were succesfully read.
    fclose(fastaFile);
    return true;
}

bool
FileLoader::loadSFFFile(const std::string& toolName, const char *fileName,
                        const bool maskBases, const bool randomizeNbases) {
    if (FileLoader::IgnoreFileName == fileName) {
        // The user does not want to use a file name. This possibly
        // happens when using peace libary from custom applications.
        return true;
    }
    SFFReader sff(fileName);
    if (!sff.isValid()) {
        std::cerr << toolName << "(Rank: ";
        std::cerr << MPI_GET_RANK()
                  << "): Error opening SFF file "
                  << fileName << " for reading." << std::endl;
        return false;
    }
    // Repeatedly read EST's from the file.
    while (sff.getPendingReads() > 0) {
        EST *est = sff.getNextRead(maskBases, randomizeNbases);
        if ((est == NULL) || (!sff.isValid())) {
            // An error occured when reading EST.
            std::cerr << toolName << ": Error loading reads from "
                      << "SFF file "  << fileName << std::endl;
            return false;
        }
    }
    // Done processing the SFF file.
    return true;
}

#endif
