#ifndef SFF_READER_CPP
#define SFF_READER_CPP

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

#include "SFFReader.h"
#include <iostream>

// A helper macro to call the combine method to combine an array of
// two bytes (unsigned char) into a short integer assuming that the
// data is in big endian order.
#define TO_SHORT(bigEndianBytePtr)                                      \
    combine<unsigned short, unsigned char>(bigEndianBytePtr[0],         \
                                           bigEndianBytePtr[1])

// A helper macro to call the combine method to combine an array of
// four bytes (unsigned char) into an unsigned integer assuming that
// the data is in big endian order.
#define TO_INT(bigEndianBytePtr)                                        \
    combine<unsigned int, unsigned short>(TO_SHORT(bigEndianBytePtr),   \
                                          TO_SHORT((bigEndianBytePtr + 2)))

// A helper macro to call the combine method to combine an array of
// eight bytes (unsigned char) into an unsigned long long assuming
// that the data is in big endian order.
#define TO_LONG_LONG(bigEndianBytePtr)                                  \
    combine<unsigned long long, unsigned int>(TO_INT(bigEndianBytePtr), \
                                              TO_INT((bigEndianBytePtr + 4)))

// Helper macro to compute the padding bytes needed so that the number
// will be an integral multiple of 8.  The last modulo with 8 handles
// the case where x is already a multiple of 8 to get the padding
// bytes to be zero.
#define PAD(x) ((8 - ((x) % 8)) % 8)

SFFReader::SFFReader(const std::string& fileName) :
    InputFile(fileName),
    sffFile(fileName.c_str(), std::ifstream::in | std::ifstream::binary) {
    // Initialize instance variables to default values
    hasValidSFFHeader = false;
    pendingReads      = 0;
    readCount         = 0;
    flowsPerRead      = 0;
    indexOffset       = 0;
    indexSize         = 0;
    // If our input stream is valid try and read the SFF header.
    if (sffFile.good()) {
        loadSFFHeader();
    }
}

SFFReader::~SFFReader() {
    sffFile.close();
}

void
SFFReader::loadSFFHeader() {
    // Flag the file as being invalid to make logic and immediate exit
    // easier in the code below.
    hasValidSFFHeader = false;
    // We assume that the file pointer is at the beginning of the
    // file.  If not, the following seekg call would be needed.
    sffFile.seekg(0, std::ifstream::beg);
    // Read the standard 31-bytes that the SFF file absolutely must
    // have as a part of the SFF header.
    unsigned char rawHeader[31];
    sffFile.read(reinterpret_cast<char*>(rawHeader), sizeof(rawHeader));
    if (sffFile.eof() || !sffFile.good()) {
        // The file does not have sufficient information to be valid.
        return;
    }
    // Validate the contents of the header to ensure it is as
    // expected.  The header must start with the magic_number field
    // value is 0x2E736666, the uint32_t encoding of the string ".sff"
    const unsigned int signature = TO_INT(rawHeader);
    if (signature != 0x2E736666) {
        // Does not have valid .sff signature.
        return;
    }
    // Verify version is 1.
    const int version = TO_INT((rawHeader + 4));
    if (version != 1) {
        // Does not have valid version.
        return;
    }
    // Save the index offset and size so that we can skip over the
    // index when loading sequences
    indexOffset = TO_LONG_LONG((rawHeader + 8));
    indexSize = TO_INT((rawHeader + 16));
    // Save the number of reads for future use.
    readCount    = TO_INT((rawHeader + 20));
    pendingReads = readCount;
    // Save the number of flows per read for future use.
    flowsPerRead = TO_SHORT((rawHeader + 28));
    // Now skip the unread bytes in the header and position the file
    // pointer at the first read header section.
    const unsigned short headerSize = TO_SHORT((rawHeader + 24));
    sffFile.seekg(headerSize, std::ifstream::beg);
    if (sffFile.eof() || !sffFile.good()) {
        // The file does not have sufficient information to be valid.
        return;
    }    
    // The SFF file seems valid.
    hasValidSFFHeader = true;
}

bool
SFFReader::isSFF(const std::string& fileName) {
    // Create a temporary reader object to make life easier.
    SFFReader tempReader(fileName);
    return tempReader.good();
}

// A macro to streamline the code in the method below
#define VALIDATE_SFF()                      \
    if (sffFile.eof() || !sffFile.good()) { \
        pendingReads = 0;                   \
        hasValidSFFHeader = false;          \
        return false;                       \
    }

bool
SFFReader::readEntry(std::string& info, std::string& sequence,
                     QualityVector& quality) {
    if (!good() || (pendingReads <= 0)) {
        // Can't read more data from the SFF file.
        return false;
    }
    // Skip over the index if we are at the index position.
    if ((unsigned) sffFile.tellg() == indexOffset) {
        // Skip over the index.
        sffFile.seekg(indexSize, std::ifstream::cur);
    }
    // First read the fixed size portion of the read header.
    unsigned char rawReadHeader[16]; 
    sffFile.read(reinterpret_cast<char*>(rawReadHeader), sizeof(rawReadHeader));
    VALIDATE_SFF();
    // Validate name length and base pair sequence lengths.
    const unsigned short nameLen = TO_SHORT((rawReadHeader + 2));
    const unsigned int   seqLen  = TO_INT((rawReadHeader + 4));
    if ((nameLen > 16000) || (seqLen > 32000)) {
        // Something is fishy here. Name appears to be longer than 16K!
        sffFile.clear(std::ifstream::failbit | std::ifstream::badbit);
        VALIDATE_SFF();
    }
    // Now extract the name for this read
    info.resize(nameLen + 1, 0);
    sffFile.read(&info[0], nameLen);
    VALIDATE_SFF();
    info[nameLen] = '\0'; // Ensure it is proper C-string
    
    // Move the file pointer to the beginning of the read data section
    // by skipping over any padding bytes that may be present. In
    // addition, we skip over the flowgram values. The number of bytes
    // of flowgram values = sizeof(unsigned short) * flowsPerRead and
    // we need to skip the flow_index_per_base array of size seqLen
    const int padSize = PAD((16 + nameLen));
    sffFile.seekg(padSize + (2 * flowsPerRead) + seqLen, std::ifstream::cur);
    VALIDATE_SFF();
    // Read the actual sequence.
    sequence.resize(seqLen + 1);
    sffFile.read(&(sequence[0]), seqLen);
    VALIDATE_SFF();
    sequence[seqLen] = '\0'; // Ensure it is proper C-string
    // Next read the quality scores into the vecotr.
    quality.resize(seqLen, 0);
    sffFile.read(reinterpret_cast<char*>(&(quality[0])), seqLen);
    // Skip the byte padding at the end. For this first we calculate
    // the padding so that it aligns to 8-byte boundary.
    const int dataPad = PAD((2 * flowsPerRead) + (3 * seqLen));
    sffFile.seekg(dataPad, std::ifstream::cur);
    VALIDATE_SFF();
    // Track the number of sequences already read
    pendingReads--;

    // Now clip out adapter sequences and low quality reads from the
    // sequence using the data in the rawReadHeader:
    const int clipQualLeft     = TO_SHORT((rawReadHeader + 8));
    const int clipQualRight    = TO_SHORT((rawReadHeader + 10));
    const int clipAdapterLeft  = TO_SHORT((rawReadHeader + 12));
    const int clipAdapterRight = TO_SHORT((rawReadHeader + 14));
    // The rules for handling the above quality and adapter ranges are
    // a bit invovled but documented at:
    // http://www.ncbi.nlm.nih.gov/Traces/trace.cgi?cmd=show&f=formats&m=doc&s=formats#sff
    const int leftIndex  = std::max(1, std::max(clipQualLeft,
                                                clipAdapterLeft)) - 1;
    const int rightIndex = std::min((clipQualRight == 0 ?
                                     seqLen : clipQualRight),
                                    (clipAdapterRight == 0 ?
                                     seqLen : clipAdapterRight));
    // Now only retain characters and quality values between leftIndex
    // and rightIndex
    sequence = sequence.substr(leftIndex, rightIndex + 1);
    sequence.reserve(rightIndex + 1); // Ensure we have space for '\0'
    sequence[rightIndex] = '\0';      // Logically terminate C-string    
    quality  = QualityVector(quality.begin() + leftIndex,
                             quality.begin() + rightIndex);
    // Everything went well.
    return true;
}

size_t
SFFReader::getCurrPos() {
    return sffFile.tellg();
}

bool
SFFReader::setCurrPos(const size_t position) {
    sffFile.seekg(position, std::ios::beg);
    return good();
}

#endif
