#ifndef BATON_ALIGNMENT_INFO_CPP
#define BATON_ALIGNMENT_INFO_CPP

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

#include "BatonAlignmentInfo.h"
#include "Utilities.h"

#include <sstream>

void
BatonAlignmentInfo::setInfo(const int estIdx, const int refESTIdx,
                            const int alignScore, const bool revCompFlag) {
    estIndex    = estIdx;
    refESTindex = refESTIdx;
    score       = alignScore;
    rcAlignment = revCompFlag;
}

void
BatonAlignmentInfo::adjustRefOffsets(const int delta) {
    for(size_t i = 0; (i < segmentList.size()); i++) {
        segmentList[i].adjustRefOffset(delta);
    }
}

void
BatonAlignmentInfo::prettyPrint(std::ostream& os, const std::string& seq,
                           const int refESTpos) const {
    int prevEndPos = -1;
    for(size_t i = 0; (i < segmentList.size()); i++) {
        const Segment& seg = segmentList[i];
        const int startPos = seg.getRefESTOffset() + refESTpos;
        ASSERT( startPos > prevEndPos );
        if (startPos > (prevEndPos + 1)) {
            // Print some blank spaces to ensure sequence is printed at the
            // correct column.
            os << std::string(startPos - prevEndPos - 1, ' ');
        }
        os << seq.substr(seg.getOthESTOffset(), seg.getLength());
        prevEndPos = startPos + seg.getLength() - 1;
    }
}

std::string
BatonAlignmentInfo::getCIGAR(const std::string& seq) const {
    ASSERT ( !segmentList.empty() );
    // End position in reference/contig    
    int prevRefEndPos = segmentList[0].getRefESTOffset();
    // End position in cDNA/EST being aligned
    int prevDnaEndPos = 0;
    // Iterate over each segment and generate CIGAR string.
    std::ostringstream cigar;
    for(size_t i = 0; (i < segmentList.size()); i++) {
        const Segment& seg = segmentList[i];
        // Compute difference between previous-end positions and
        // current-start positions to streamline detection of indels
        // etc. (right below in series of if-statements)
        const int refDelta = seg.getRefESTOffset() - prevRefEndPos;
        const int dnaDelta = seg.getOthESTOffset() - prevDnaEndPos;
        ASSERT ( refDelta >= 0 );  // Deltas must never be negative
        ASSERT ( dnaDelta >= 0 );  // in any situation.
        if (prevDnaEndPos > 0) {
            // Use difference in previous-end position and current-start
            // positions to detect type of operation.
            if ((dnaDelta == 0) && (refDelta > 0)) {
                // The segments in other-cDNA are consecutive without
                // gaps. So it must be an delete.
                cigar << refDelta << "D";
            } else if ((refDelta == 0) && (dnaDelta > 0)) {
                // The segments in contig are consecutive without gaps. So
                // it must be an insert.
                cigar << dnaDelta << "I";
            } else if ((refDelta > 0) && (dnaDelta > 0)) {
                // There is a region of mismatch between the two. We
                // could try and detect additional
                // inserts/delets; but for now we just report mismatch
                cigar << dnaDelta << "X";
            }
        } else {
            // This is the first segment. No previous segments are
            // available for comparisons. Generate initial clippings
            // if any.
            if (seg.getOthESTOffset() > 0) {
                // The first several nucleotides are being soft-clipped
                cigar << seg.getOthESTOffset() << "S";
            }
        }
        // Generate the CIGAR sub-string for matching nucleotides
        cigar << seg.getLength() << "M";
        // Update previous-end position trackers
        prevRefEndPos = seg.getRefESTEndOffset() + 1;
        prevDnaEndPos = seg.getOthESTEndOffset() + 1;
    }
    // Finally create a hard clip for trailing nucleotide sequences as
    // needed.    
    const int rightOverhang = seq.size() - prevDnaEndPos;
    if (rightOverhang > 0) {
        // We have some overhang that are soft-clipped
        cigar << rightOverhang << "S";
    }
    // Return the CIGAR string back to the caller
    return cigar.str();
}

#endif
