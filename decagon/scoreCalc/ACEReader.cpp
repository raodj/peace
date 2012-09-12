#ifndef ACE_READER_CPP
#define ACE_READER_CPP

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

#include "ACEReader.h"
#include "Utilities.h"
#include "ESTList.h"
#include "Contig.h"

#include <sstream>

ACEReader::ACEReader(const std::string& fileName) throw(DecagonException) {
    // Open ACE file.
    aceFile.open(fileName.c_str());
    if (!aceFile.good()) {
        throw EXP(100, "Error opening supplied file",
                  "Ensure file name and path is valid."); 
    }
    // If open is successful, try and read the AS header
    std::string asHeader;
    getNextLine(asHeader);
    // Check to ensure we have a valid AS header
    if (!aceFile.good() || (asHeader.find("AS") != 0)) {
        throw EXP(101, "Unable to read ACE header (AS) from file.",
                  "Ensure file is actually an ACE file.");
    }
}

ACEReader::~ACEReader() {
    aceFile.close();
}

bool
ACEReader::getNextLine(std::string& line) {
    if (!aceFile.good()) {
        return false;
    }
    // Clear out existing line contents if any
    line.clear();
    do {
        // Read next-line from file.
        if (std::getline(aceFile, line).good()) {
            // Trim-off any leading and trailing white spaces.
            trim(line);
        }
    } while ((aceFile.good()) && (line.empty()));
    // Read a line?
    return !line.empty();
}

bool
ACEReader::readContig(Contig& contig, ESTList& reads)
    throw (DecagonException) {
    // Read contig CO header and and set the contig ID
    const std::string contigID = readCOHeader();
    if (contigID == "") {
        // End of file reached! Nothing more to do
        return false;
    }
    contig.setId(contigID);
    // Read the contig consensus sequence
    std::string consensus;
    readMultiLineData(consensus);
    contig.setConsensus(consensus);
    // Create empty nucleotide occurrence frequency vector for this
    // contig.
    NtDistrList ntDistr(consensus.size());
    // Check and read the BQ header and set quality values
    checkReadContigQuality(contig);
    // Next load AF and BS entries for this contig.
    StringIntMap afEntries;
    BSLookupMap  bsEntries;
    checkReadAFBSEntries(afEntries, bsEntries);
    // Next load all the reads for this contig from the file.
    while (loadRead(contig, reads, ntDistr, afEntries));
    // Next read and skip all RT, CT, and WA entries for this contig.
    while (readCtRtWaEntries());
    // Finally, load any trailing reads for this contig from the file.
    // Mira generates reads aft RT and CT entries even though ACE file
    // format stipulates reads should appear before RT, CT, and WA
    // entries. So this is a work around for a Mira bug.
    while (loadRead(contig, reads, ntDistr, afEntries));    
    // Contig read successfully
    return true;
}

bool
ACEReader::readCtRtWaEntries() throw (DecagonException) {
    // Check to see if the next line has CT or RD entry. If so, read it.
    std::ifstream::streampos savedFilePos = aceFile.tellg();
    std::string headerLine;
    getNextLine(headerLine);
    const std::string header = headerLine.substr(0, 2);
    if ((header != "RT") && (header != "CT") && (header != "WA")) {
        // This is not an RD or CT entry. Don't process it.
        aceFile.seekg(savedFilePos);
        return false;
    }
    // This is a CT/RT/WA entry. These headers can have nested "{" in
    // some cases. The header line has the first "{".
    ASSERT ( headerLine.find('{') != std::string::npos );
    int backCount = 1;
    while (backCount > 0) {
        if (!getNextLine(headerLine)) {
            // Error occured when reading header.
            throw EXP(112, "Error when reading RT/CT/WA entries",
                      "The ACE file is not valid");
        }
        if (headerLine.find('{') != std::string::npos) {
            backCount++;
        } else if (headerLine.find('}') != std::string::npos) {
            backCount--;
        }
    }
    // Successfully read a CT/RT/WA entry.
    return true;
}

bool
ACEReader::loadRead(Contig& contig, ESTList& reads, NtDistrList& ntDistr,
                    StringIntMap& afEntries)
    throw (DecagonException) {
    // Check to see if the next line has RD entry. If so, read it.
    std::ifstream::streampos savedFilePos = aceFile.tellg();
    std::string rdHeader;
    getNextLine(rdHeader);
    if (rdHeader.substr(0, 3) != "RD ") {
        // This is not an RD entry. Don't process it.
        aceFile.seekg(savedFilePos);
        return false;
    }
    // This is an RD entry. Obtain the read name from the header.
    // Extract the contig ID from the header.
    const size_t nonSpacePos = rdHeader.find_first_not_of(' ', 2);
    const size_t spcAfterID  = rdHeader.find(' ', nonSpacePos);
    const std::string readID = rdHeader.substr(nonSpacePos,
                                               spcAfterID - nonSpacePos);
    // Read the nucleotide sequence for the entry.
    std::string nucleotides;
    readMultiLineData(nucleotides);
    // Add nucleotide to the read list.
    const int readIndex      = reads.size();
    reads.add(readIndex, readID, nucleotides);
    // Tally the nucleotide occurrences in the contig.
    StringIntMap::const_iterator posValue = afEntries.find(readID);
    if (posValue == afEntries.end()) {
        // Did not find the AF entry for this read!
        std::string errMsg = "Did not find AF entry for read: ";
        errMsg += readID;
        throw EXP(111, errMsg.c_str(), "The ACE file is not valid");
    }
    // In ACE file index position start with 1. So we need -1 below.
    const int startPos = posValue->second - 1;
    ASSERT ( startPos >= 0 );
    for(size_t i = 0; (i < nucleotides.size()); i++) {
        ntDistr[i + startPos].add(nucleotides[i]);
    }
    
    // Read QA header if present.
    savedFilePos = aceFile.tellg();
    std::string qaHeader;
    getNextLine(qaHeader);
    if (qaHeader.substr(0, 3) != "QA ") {
        // This is not a QA line. Undo line reading
        aceFile.seekg(savedFilePos);
    }
    // Read DS header if present.
    savedFilePos = aceFile.tellg();
    std::string dsHeader;
    getNextLine(dsHeader);
    if (dsHeader.substr(0, 3) != "DS ") {
        // This is a DS line. Undo line reading
        aceFile.seekg(savedFilePos);
    }

    // Finally create and add alignment information.
    AlignmentInfo ai(readIndex, false, startPos, "", 255, 0, 0);
    contig.addAlignmentInfo(ai);
    // Successfully read an RD entry.
    return true;
}

void
ACEReader::checkReadAFBSEntries(StringIntMap& afEntries,
                                BSLookupMap& bsEntries)
    throw (DecagonException) {
    // Read and process entires until we encounter some other entry.
    std::string line;
    do {
        const std::ifstream::streampos savedFilePos = aceFile.tellg();
        if (!getNextLine(line)) {
            // Error reading next line.
            break;
        }
        // Wrap line into a stream to ease reading of fields
        std::istringstream infoStream(line);
        std::string header;
        infoStream >> header; // Extract the header out for convenience.
        if (header == "AF") {
            // Extract Read-ID, U/C flag, and start position
            std::string readID, compFlag;
            int startPos = -1;
            infoStream >> readID >> compFlag >> startPos;
            ASSERT(startPos > 0);
            // Add information to the afEntries hash map.
            afEntries[readID] = startPos;
        } else if (header == "BS") {
            // Extract startPos, endPos, readID.
            int startPos, endPos;
            std::string readID;
            infoStream >> startPos >> endPos >> readID;
            // Add information to the bsEntries hash map.
            bsEntries[readID] = std::make_pair(startPos, endPos);
        } else if (!line.empty()) {
            // This is not an AF or a BS entry. Reset file pointer and
            // break out.
            aceFile.seekg(savedFilePos);
            break;
        }
        // This do-while loop breaks when we encounter a non-blank
        // line that is not an AF or BS line.
    } while (aceFile.good());
}

void
ACEReader::checkReadContigQuality(Contig& contig)
    throw (DecagonException) {
    // Contig quality may not be present. So here we tentatively read
    // the next line looking for "BQ" entry. If it is not found then
    // we reset file pointer back and return.
    const std::ifstream::streampos savedFilePos = aceFile.tellg();
    // Read next line and check if it is "BQ". If not reset file
    // pointer and return.
    std::string bqHeader;
    getNextLine(bqHeader);
    if (bqHeader != "BQ") {
        aceFile.seekg(savedFilePos);
        return;
    }
    // We found the BQ header. Read multiple-lines (until a blank-line
    // is read) to obtain quality values.
    std::string qualString;
    readMultiLineData(qualString, " ");
    // Convert quality string into integers to build quality
    // vector. For "*" consensus characters, there will not be a
    // quality read and the entries need to be appropirately handled
    // in the loop below.
    const std::string consensus = contig.getConsensus();
    std::istringstream qualStream(qualString);
    std::vector<int> qualList;
    for(size_t i = 0; ((i < consensus.size()) && !qualStream.eof()); i++) {
        if (consensus[i] == '*') {
            // There is no quality read for '*'
            qualList.push_back(-1);
        } else {
            // Read next base quality value and use it.
            int baseQuality;
            qualStream >> baseQuality;
            qualList.push_back(baseQuality);
        }
    }
    // Ensure we have the correct number of values (sanity check)
    if (contig.getConsensus().size() != qualList.size()) {
        std::ostringstream errMsg;
        errMsg << "The number of nucleotides in the consensus did not "
               << "match number of quality values read [consensus length="
               << contig.getConsensus().size() << ", base quality count="
               << qualList.size() << "]";
        throw EXP(106, errMsg.str().c_str(),
                  "Either the ACE file is bad or it is incompatible");
    }
    // Set the quality values for the contig.
    contig.setQuality(qualList);
}

bool
ACEReader::readMultiLineData(std::string& data, const std::string& eol) {
    // Get the next non-blank line.
    if (!getNextLine(data)) {
        // Can't even get one non-blank line.
        return false;
    }
    // Next, repeatedly read lines until we get a blank-line.
    std::string nextLine;
    do {
        // Read next-line from file.
        if (std::getline(aceFile, nextLine).good()) {
            // Trim-off any leading and trailing white spaces.
            trim(nextLine);
            // Concatnate the results.
            data += eol;
            data += nextLine;
        }
    } while ((aceFile.good()) && (!nextLine.empty()));
    // Read a line?
    return !data.empty();
}

std::string
ACEReader::readCOHeader() throw(DecagonException) {
    // Read the CO line from the file.
    std::string coHeader;
    getNextLine(coHeader);
    if ((aceFile.eof()) && (coHeader.empty())) {
        // No more contigs to process
        return "";
    }
    if (!aceFile.good() || (coHeader.find("CO") != 0)) {
        throw EXP(102, "Unable to read contig header (CO) from file.",
                  "Ensure file is actually an ACE file.");
    }
    // Extract the contig ID from the header.
    const size_t nonSpacePos = coHeader.find_first_not_of(' ', 2);
    const size_t spcAfterID  = coHeader.find(' ', nonSpacePos);
    std::string contigID     = coHeader.substr(nonSpacePos,
                                               spcAfterID - nonSpacePos);
    // Return contigID back to the caller (we ignore the other fields
    // for now).
    return contigID;
}

void
ACEReader::testReader(const std::string& aceFileName) {
    try {
        ACEReader ar(aceFileName);
        Contig contig;
        ESTList reads;
        while (ar.readContig(contig, reads));
    } catch (const DecagonException& de) {
        std::cerr << de << std::endl;
    }
}

#endif
