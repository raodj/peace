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
#include "ESTCodec.h"
#include <string>
#include <algorithm>
#include <iostream>
#include <cctype>
#include <cstring>
#include <cstdio>
#include <cstdlib>

EST::EST(const int idValue, const std::string& information,
         const std::string& seq) :
    id(idValue), info(information), sequence(seq),
    sequenceLen((int) seq.size()),  customData(NULL) {
    similarity = 0;
    processed  = false;
    populated  = true;
}

EST::EST(const int idValue, const std::string& information,
         const std::string& seq, const QualityVector& quals) :
    id(idValue), info(information), sequence(seq), 
    sequenceLen((int) seq.size()), quality(quals), customData(NULL) {
    similarity = 0;
    processed  = false;
    populated  = true;
}

EST::EST(const EST& src) :
    id(src.id), sequenceLen(src.sequenceLen) {
    info       = src.info;
    sequence   = src.sequence;
    quality    = src.quality;
    similarity = src.similarity;
    processed  = src.processed;
    customData.reset(src.customData.get());
    populated  = src.populated;
}

EST::~EST() {
    unpopulate();
}

void
EST::unpopulate() {
    // Clear out the information and reduce memory allocated to zero.
    info.resize(0);
    sequence.resize(0);
    quality.resize(0);
    // Set populated flag to false to indicate that we don't have all
    // the data for this EST.
    populated = false;
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

void
EST::repopulate(const std::string& info, const std::string& sequence,
                const QualityVector& quality) {
    this->info      = info;
    this->sequence  = sequence;
    this->quality   = quality;
    this->populated = true;
}

// A dummy default constructor to keep the compiler happy.
EST::EST() : id(-1), info(NULL), sequence(NULL), sequenceLen(0), 
             similarity(0), populated(false) {}

// A dummy operator=
EST&
EST::operator=(const EST&) {
    return *this;
}

std::string
EST::normalizeBases(const std::string& srcSequence,
                    const bool maskBases, const bool randomizeNbases,
                    const int randSeed) {
    std::string sequence(srcSequence);
    const size_t seqLen = sequence.size();
    static const std::string LowCaseBases = "atcg";
    static const std::string UpCaseBases  = "ATCG";
    // Setup random number generator seed for use further below.
    if (randomizeNbases) {
        srandom(randSeed);
    }
    // Normalize the sequence.
    for(size_t i = 0; (i < seqLen); i++) {
        size_t index = 0;
        // Obtain base to be normalized.
        char nt   = sequence[i];
        if ((index = LowCaseBases.find(nt)) != std::string::npos) {
            nt = (maskBases ? 'N' : UpCaseBases[index]);
        } else if (UpCaseBases.find(nt) == std::string::npos) {
            // An entry other than "atcgATCG" to be ignored
            nt = 'N';
        }
        // Check and randomize 'N' bases as directed.
        if (randomizeNbases && (nt == 'N')) {
            // Randomly select a base here.
            register int ntIndex = (random() >> 5) & 3;
            nt = UpCaseBases[ntIndex];
        }
        // Store normalize character
        sequence[i] = nt;
    }
    // Return the normalized sequence back
    return sequence;
}

std::string
EST::getRCSequence() const {
    // Create a reversed-copy of the incoming string.
    std::string rc(sequence.rbegin(), sequence.rend());
    // Flip each nucleotide in the sequence to its complement.
    for(size_t idx = 0; (idx < rc.size()); idx++) {
        rc[idx] = ESTCodec::getComp(rc[idx]);
    }
    // Return the reverse complement version.
    return rc;
}

#endif
