#ifndef FASTA_FILE_CPP
#define FASTA_FILE_CPP

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

#include "FASTAFile.h"

FASTAFile::FASTAFile(const std::string& fileName) : InputFile(fileName) {
    // Open the FASTA file to be read and processed
    fastaFile.open(fileName.c_str());
}

FASTAFile::~FASTAFile() {
    fastaFile.close();
}

bool
FASTAFile::isFASTA(const std::string& fileName) {
    // Create a temporary FASTAFile object to make life easier
    FASTAFile tempFile(fileName);
    // Read the first entry to validate contents.
    std::string info, sequence;
    QualityVector quality;
    const bool result = tempFile.readEntry(info, sequence, quality);
    // The result of trying to read a sequence provide sufficiently
    // reasonable confidence that this file is a FASTA file.
    return result;
}

bool
FASTAFile::hasNextEntry() {
    return !fastaFile.eof();
}

size_t
FASTAFile::getCurrPos() {
    return fastaFile.tellg();
}

bool
FASTAFile::setCurrPos(const size_t position) {
    fastaFile.seekg(position, std::ios::beg);
    if (!fastaFile.fail()) {
        // The seek operation was successful. If the stream was not in
        // a good state earlier (because we hit EOF), we reset the EOF
        // flag here to permit smooth reading after this seek.
        fastaFile.clear();
        return true;
    }
    // Error during seek
    return false;
}

bool
FASTAFile::readEntry(std::string& info, std::string& sequence,
                     QualityVector& quality) {
    // Some basic validation on the fastFile first.    
    if (!fastaFile.good()) {
        // Can't read information from the fasta file.
        return false;
    }
    // Check if the current character in the fasta file is a ">"
    if (fastaFile.peek() != '>') {
        // No! This is not a valid fasta file header.
        return false;
    }
    // Ok, now read rest of the line as the information from the file.
    // The information should not exceed 1024 characters as per the
    // standard?  Read header line into info.
    std::getline(fastaFile, info);
    // Now read the actual sequence information from the fasta file.
    // This is performed until either EOF is reached or the next
    // header is detected.
    sequence.clear();
    int headerChar;    
    do {
        // Check if the next char is a header start without consuming
        // it permanently.
        if ((headerChar = fastaFile.peek()) != EOF) {
            if (headerChar != '>') {
                // This is still sequence information. Read it in.
                std::string ntSeq;
                std::getline(fastaFile, ntSeq);
                sequence += ntSeq;
            }
        }
    } while (fastaFile.good() && (headerChar != '>'));
    // Simply clear out the quality value as the FASTA file does not
    // provide any.
    quality.clear();
    // Return success as everything went well.
    return !fastaFile.fail();
}

#endif
