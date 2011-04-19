#ifndef ALIGNMENT_ALGORITHM_CPP
#define ALIGNMENT_ALGORITHM_CPP

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
// Authors:   Jenna Zhang               zhangy9@muohio.edu
//            Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include "AlignmentAlgorithm.h"
#include "Utilities.h"
#include <cstring>
#include <algorithm>

AlignmentAlgorithm::AlignmentAlgorithm(const int match,
                                       const int mismatch,
                                       const int gap,
                                       const int aScoreDelPenParm,
                                       const int aScoreInsPenParm,
                                       const int aScoreSubPenParm)  {
    // We create a 2-D matrix of scores to make alignment faster. The
    // matrix is initialized by helper method.
    computeScoreMatrix(match, mismatch);
    gapPenalty = gap;
    // Setup the a-score parameters
    aScoreDelPen = aScoreDelPenParm;
    aScoreInsPen = aScoreInsPenParm;
    aScoreSubPen = aScoreSubPenParm;
    // Setup the base encodings using the order of bases in the string.
    const std::string BaseOrder = "ACGTNPYRWSKMDVHBX";
    for(size_t idx = 0; (idx < BaseOrder.size()); idx++) {
        baseCode[(int) BaseOrder[idx]] = idx;
    }
}

void
AlignmentAlgorithm::computeScoreMatrix(const int match, const int mismatch) {
    // We create a 2-D matrix of scores to make alignment faster. The
    // matrix is initialized below in a two-step process to provide
    // maximum compatibility between C++ compilers.  First we
    // introduce m (match) and n (mismatch) homonyms to make the
    // matrix look nice in code.
	const int m = match;
	const int n = mismatch;
    
	const int tmpScoringMatrix[17][17] = {
        //A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X        
        {m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n},  //A
        {n, m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n},  //C
        {n, n, m, n, n, n, n, n, n, n, n, n, n, n, n, n, n},  //G
        {n, n, n, m, n, n, n, n, n, n, n, n, n, n, n, n, n},  //T
        {n, n, n, n, m, n, n, n, n, n, n, n, n, n, n, n, n},  //N
        {n, n, n, n, n, m, n, n, n, n, n, n, n, n, n, n, n},  //P: gap
        {n, m, n, m, n, n, m, n, n, n, n, n, n, n, n, n, n},  //Y: C, T
        {m, n, m, n, n, n, n, m, n, n, n, n, n, n, n, n, n},  //R: A, G
        {m, n, n, m, n, n, n, n, m, n, n, n, n, n, n, n, n},  //W: A, T
        {n, m, m, n, n, n, n, n, n, m, n, n, n, n, n, n, n},  //S: G, C
        {n, n, m, m, n, n, n, n, n, n, m, n, n, n, n, n, n},  //K: T, G
        {m, m, n, n, n, n, n, n, n, n, n, m, n, n, n, n, n},  //M: C, A
        {m, n, m, m, n, n, n, n, n, n, n, n, m, n, n, n, n},  //D: not C
        {m, m, m, n, n, n, n, n, n, n, n, n, n, m, n, n, n},  //V: not T
        {m, m, n, m, n, n, n, n, n, n, n, n, n, n, m, n, n},  //H: not G
        {n, m, m, m, n, n, n, n, n, n, n, n, n, n, n, m, n},  //B: not A
        {n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, m}   //X: unknown
    };

    // Copy the data from local variables to our instance variables.
    ASSERT( sizeof(scoreMatrix) == sizeof(tmpScoringMatrix) );
    memcpy(scoreMatrix, tmpScoringMatrix, sizeof(scoreMatrix));
}

AlignmentAlgorithm::~AlignmentAlgorithm() {
    // Currently, this class does not use any dynamic memory. So the
    // destructor is empty.
}

void
AlignmentAlgorithm::setAlignmentParams(const int match, const int mismatch,
                                       const int gap) {
    computeScoreMatrix(match, mismatch);
    gapPenalty = gap;    
}

void
AlignmentAlgorithm::setAScoreParams(const int aScoreDelPen,
                                    const int aScoreInsPen,
                                    const int aScoreSubPen) {
    this->aScoreDelPen = aScoreDelPen;
    this->aScoreInsPen = aScoreInsPen;
    this->aScoreSubPen = aScoreSubPen;    
}

int
AlignmentAlgorithm::getNWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   int& alignScore,
                                   int& aScore,
                                   std::string& alignedSeq1,
                                   std::string& alignedSeq2) const {
	//calculate the alignment score and record the traceback path
	int numOfRows = seq1.length();
	int numOfCols = seq2.length();
    // Setup the alignment and trace matrices
    int alignMatrix[numOfRows + 1][numOfCols + 1];
	DIR trace[numOfRows + 1][numOfCols + 1];
	trace[0][0]       = BAD_DIR;
    alignMatrix[0][0] = 0;
	//initialize the first row and column of the score and trace
	//matrices as per needleman-wunch algorithm
	for (int col = 1; col <= numOfCols; col++) {
		alignMatrix[0][col] = gapPenalty * col;
		trace[0][col]       = WEST;
	}
	for (int row = 1; row <= numOfRows; row++) {
		alignMatrix[row][0] = gapPenalty * row;
		trace[row][0]       = NORTH;
	}

	// Build the needleman-wunch dynamic programming matrix row-by-row
	for (int row = 1; row <= numOfRows; row++) {
		for (int col = 1; col <= numOfCols; col++) {
			// Initialize max to the first of the three terms (NORTH).
			DIR dirFlag  = NORTH;
			int maxScore = alignMatrix[row - 1][col] + gapPenalty;

			// See if the second term is larger (WEST).
			const int westScore = alignMatrix[row][col - 1] + gapPenalty;
			if (maxScore <= westScore) {
				maxScore = westScore;
				dirFlag  = WEST;
			}

			// See if the third term is the largest (NORTHWEST)
			const int base1 = encodeBase(seq1[row - 1]);  // Short cut
			const int base2 = encodeBase(seq2[col - 1]);  // Short cut
			const int northWestScore = scoreMatrix[base1][base2] +
                alignMatrix[row - 1][col - 1];
			if (maxScore <= northWestScore) {
				maxScore = northWestScore;
				dirFlag  = NORTH_WEST;
			}
            // Save the best score and direction of move in matrix for
            // use to reconstruct the aligned sequences below.
			alignMatrix[row][col] = maxScore;
			trace[row][col]       = dirFlag;
		}
	}
    // Now store the final alignment score in the parameter
	alignScore = alignMatrix[numOfRows][numOfCols];

	// Trace back and get the alignment strings. Note that the
	// nucleotide sequence in tStr1 and tStr2 are reverse as we build
	// the alignments below.
    std::string tStr1;
    std::string tStr2;
    // Reserve sufficient memory in tStr1 and tStr2 to improve efficiency.
    const size_t maxLen = std::max(seq1.size(), seq2.size());
    tStr1.reserve(maxLen);
    tStr2.reserve(maxLen);
    
	int row = numOfRows;
	int col = numOfCols;
    // Compute a-score as we go along constructing the aligned
    // sequences assuming seq1 is the reference sequence.
    aScore = 2 * seq1.size();
	while ((row > 0) || (col > 0)) {
		switch (trace[row][col]) {
		case NORTH: //row-1, col, north
            row--;  // Must be done first as it affects code below
			tStr1.append(1, seq1[row]);
			tStr2.append(1, '-');
            aScore += aScoreDelPen; // Deletion detected
			break;
		case WEST: //row, col-1, west
            col--; // Must be done first as it affects code below
			tStr1.append(1, '-');
			tStr2.append(1, seq2[col]);
            aScore += aScoreInsPen; // Insertion detected
			break;
		case NORTH_WEST: // row-1, col-1, northwest
            row--; // These two decrements msut be done first as they
            col--; // affect seq1 and seq2 accesses right below.
			tStr1.append(1, seq1[row]);
			tStr2.append(1, seq2[col]);
            // The nucleotides could be the same or different. 
            aScore += (seq1[row] != seq2[col] ? aScoreSubPen : 0);
			break;
            
        case BAD_DIR:
        default:
            // Do nothing
            break;
		}
	}

	// Recollect that tStr1 and tStr2 are reverse of the final
	// sequence we expect. So reverse the string and copy the result
	// into the final outgoing parameters.
    alignedSeq1.resize(tStr1.size(), '-');
    std::reverse_copy(tStr1.begin(), tStr1.end(), alignedSeq1.begin());
    // Reverse copy the second sequence into the outgoing parameter
    alignedSeq2.resize(tStr2.size(), '-');
    std::reverse_copy(tStr2.begin(), tStr2.end(), alignedSeq2.begin());
    // return the score
    return alignScore;
}

#endif
