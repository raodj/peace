#ifndef EST_CPP
#define EST_CPP

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
#include <string>
#include <algorithm>
#include <iostream>
#include <cctype>
#include <cstring>
#include <cstdio>

// The shared static list of ESTs currently available.
std::vector<EST*> EST::estList;

// The length of the longest EST we have
size_t EST::maxESTlen = 0;

// A "dummy" EST for returning in place of filtered-out ESTs
EST* dummy = new EST(-1, "", "", 0);
        
EST::EST(const int idValue, const char *information, const char* seq,
         const int fileOffset) : id(idValue), offset(fileOffset),
                                 customData(NULL) {
    // Duplicate the information and sequnece data as needed.
    info       = EST::duplicate(information);
    sequence   = EST::duplicate(seq);
    similarity = 0;
    processed  = false;
}

EST::~EST() {
    unpopulate();
}

void
EST::unpopulate() {
    if (info != NULL) {
        delete [] info;
    }
    if (sequence != NULL) {
        delete [] sequence;
    }
}

std::string
EST::getLine(FILE *fastaFile) {
    char buffer[1024];
    std::string retVal;
    char *result = NULL;
    
    do {
        // Read a line from the file.
        result = fgets(buffer, 1023, fastaFile);
        const int len = (int) strlen(buffer);
        if (buffer[len - 1] == '\n') {
            // Remove trailing newline as we don't really need it.
            buffer[len - 1] = '\0';
            // Set result to NULL to indicate termination
            result = NULL;
        }
        // Add data read from file to the string.
        retVal += buffer;
    } while (result != NULL);

    // return the string loaded from the file back to the caller.
    return retVal;
}

EST*
EST::create(const int id, const char *info,
            const char* sequence, const long offset) {
    if (id != (int) estList.size()) {
        // This id is not acceptable. Sorry.
        return NULL;
    }
    // Instantiate new EST
    EST *newEST = new EST(id, info, sequence, offset);
    // Add est to end of the est list.
    estList.push_back(newEST);
    // return the newly created EST back tot he caller
    return newEST;
}

EST*
EST::create(FILE* fastaFile, int& lineNum) {
    // Some basic validation on the fastFile first.
    if (feof(fastaFile) || ferror(fastaFile)) {
        // Can't read information from the fasta file.
        return NULL;
    }
    // Save current offset in the fasta file for future reference when
    // creating an EST further down in this method.
    const long offset = ftell(fastaFile);
    
    // Check if the current character in the fasta file is a ">"
    int headerChar;
    if ((headerChar = fgetc(fastaFile)) != '>') {
        // No! This is not a valid fasta file.
        ungetc(headerChar, fastaFile);
        return NULL;
    }
    // Ok, now read rest of the line as the information from the file.
    // The information should not exceed 1024 characters as per the
    // standard?
    // Read header line into info.
    std::string headerLine = getLine(fastaFile);
    lineNum++;
    // Now read the actual sequence information from the fasta file.
    // This is performed until either EOF is reached or the next
    // header is detected.
    std::string sequence;
    do {
        // Check if the next char is a header start without consuming
        // it permanently.
        if ((headerChar = fgetc(fastaFile)) != EOF) {
            ungetc(headerChar, fastaFile);
            if (headerChar != '>') {
                // This is still sequence information. Read it in.
                sequence += getLine(fastaFile);
                // Track line numbers
                lineNum++;
            }
        }
    } while (!feof(fastaFile) && !ferror(fastaFile) && (headerChar != '>'));
    // Check if all the data was loaded successfully
    if (!ferror(fastaFile)) {
        // Convert est information to upper case
        std::transform(sequence.begin(), sequence.end(), sequence.begin(),
                       toupper);
        const char* const seqBP = sequence.c_str();
        // Compute new max EST length
        maxESTlen = std::max(maxESTlen, strlen(seqBP));
        // Create a new est with all the information.
        EST *est = new EST((int) estList.size(), headerLine.c_str(), seqBP,
                           offset);
        // Add it to the est list.
        estList.push_back(est);
        // Return newly created est back to the caller.
        return est;
    }

    // Can't create a valid EST
    return NULL;
}

void
EST::deleteAllESTs() {
    const int ESTCount = (int) estList.size();
    for(int id = 0; (id < ESTCount); id++) {
        delete estList[id];
    }
    // Clear out all entries as they have all been deleted.
    estList.clear();
}

void
EST::deleteLastESTs(const int count) {
    for(int i = 0; ((estList.size() > 0) && (i < count)); i++) {
        // Remove the last entry
        delete estList.back();
        estList.pop_back();
    }
}

void
EST::dumpESTList(std::ostream& os) {
    const int EstCount = (int) estList.size();
    for(int id = 0; (id < EstCount); id++) {
        // Dump the EST information out.
        estList[id]->dumpEST(os);
    }
}

void
EST::dumpESTList(std::ostream& os, const bool processed) {
    const int EstCount = (int) estList.size();
    for(int id = 0; (id < EstCount); id++) {
        // Dump the EST information out based on processed status.
        if (estList[id]->processed == processed) {
            // Matched condition.
            estList[id]->dumpEST(os);
        }
    }
}

void
EST::dumpEST(std::ostream& os) {
    const int LineSize = 100;
    os << ">";
    os << getInfo() << std::endl;
    // Dump out the sequence to that no sequence display is longer
    // than LineSize characters.
    const char *seq   = getSequence();
    const int  seqLen = (int) strlen(seq);
    for(int pos = 0; (pos < seqLen); pos++) {
        if ((pos > 0) && ((pos % LineSize) == 0)) {
            os << "\n";
        }
        os << seq[pos];
    }
    os << "\n";
}

char*
EST::duplicate(const char *src) {
    char *copy = NULL;
    if (src != NULL) {
        const size_t len = strlen(src) + 1;
        copy = new char[len];
#ifdef _WINDOWS
        strcpy_s(copy, len, src);
#else
        strncpy(copy, src, len);
#endif
    }
    return copy;
}

// A dummy default constructor to keep the compiler happy.
EST::EST() : id(-1), info(NULL), sequence(NULL), offset(-1), similarity(0) {}

// A dummy operator=
EST&
EST::operator=(const EST&) {
    return *this;
}

int
EST::getProcessedESTCount() {
    int count = 0;
    const int ESTCount = (int) estList.size();
    for(int id = 0; (id < ESTCount); id++) {
        if (estList[id]->processed) {
            count++;
        }
    }
    return count;
}

size_t
EST::getMaxESTLen() {
    if (maxESTlen == 0) {
        // First time we have to compute the max length as it is not
        // yet known.
        const int ESTCount = (int) estList.size();
        for(int id = 0; (id < ESTCount); id++) {
            maxESTlen = std::max(maxESTlen, strlen(estList[id]->getSequence()));
        }
    }
    // We have the maximum length EST
    return maxESTlen;
}

#endif
