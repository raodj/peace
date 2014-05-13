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
#include <iterator>
#include <limits>

// The static a-score matrix that contains fractional values to be
// used.  Note this matrix should be symmetric
const int AlignmentAlgorithm::aScoreMatrix[17][17] = {
       //A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X        
        {1, 0, 0, 0, 1, 0, 0, 2, 2, 0, 0, 2, 3, 3, 3, 0, 1},  //A
        {0, 1, 0, 0, 1, 0, 2, 0, 0, 2, 0, 2, 0, 3, 3, 3, 1},  //C
        {0, 0, 1, 0, 1, 0, 0, 2, 0, 2, 2, 0, 3, 3, 0, 3, 1},  //G
        {0, 0, 0, 1, 1, 0, 2, 0, 2, 0, 2, 0, 3, 0, 3, 3, 1},  //T
        {1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},  //N
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},  //P: gap
        {0, 2, 0, 2, 1, 0, 1, 0, 3, 3, 3, 3, 3, 3, 2, 2, 1},  //Y: C, T
        {2, 0, 2, 0, 1, 0, 0, 1, 3, 3, 3, 3, 2, 2, 3, 3, 1},  //R: A, G
        {2, 0, 0, 2, 1, 0, 3, 3, 1, 0, 3, 3, 2, 3, 2, 3, 1},  //W: A, T
        {0, 2, 2, 0, 1, 0, 3, 3, 0, 1, 3, 3, 3, 2, 3, 2, 1},  //S: G, C
        {0, 0, 2, 2, 1, 0, 3, 3, 3, 3, 1, 0, 2, 3, 3, 2, 1},  //K: T, G
        {2, 2, 0, 0, 1, 0, 3, 3, 3, 3, 0, 1, 3, 2, 2, 3, 1},  //M: C, A
        {3, 0, 3, 3, 1, 0, 3, 2, 2, 3, 2, 3, 1, 2, 2, 2, 1},  //D: not C
        {3, 3, 3, 0, 1, 0, 3, 2, 3, 2, 3, 2, 2, 1, 2, 2, 1},  //V: not T
        {3, 3, 0, 3, 1, 0, 2, 3, 2, 3, 3, 2, 2, 2, 1, 2, 1},  //H: not G
        {0, 3, 3, 3, 1, 0, 2, 3, 3, 2, 2, 3, 2, 2, 2, 1, 1},  //B: not A
        {1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}   //X: unknown
    };


template <typename T>
void writeVector(const std::vector<T>& vec) {
    std::copy(vec.begin(), vec.end(),
              std::ostream_iterator<T>(std::cout, " "));
    std::cout << std::endl;
}

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
    // Verify that the aScoreMatrix is symmetric when developer
    // assertions are enabled.
    ASSERT( verifyAScoreMatrixIsSymmetric() );
}

bool
AlignmentAlgorithm::verifyAScoreMatrixIsSymmetric() {
    const std::string BaseOrder = "ACGTNPYRWSKMDVHBX";    
    const int MatSize = 17;
    for(int row = 0; (row < MatSize); row++) {
        for(int col = row; (col < MatSize); col++) {
            if (aScoreMatrix[row][col] != aScoreMatrix[col][row]) {
                // Found a mismatch!
                std::cerr << "Warning:: Asymmetry detected: "
                          << "aScoreMatrix[" << BaseOrder[row] << ","
                          << BaseOrder[col]  << "]="
                          << aScoreMatrix[row][col] << " but "
                          << "aScoreMatrix[" << BaseOrder[col] << ","
                          << BaseOrder[row]  << "]="
                          << aScoreMatrix[col][row] << std::endl;
            }
        }
    }
    // This method always returns true
    return true;
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
                                   int& alignScore, int& aScore,
                                   std::string& alignedSeq1,
                                   std::string& alignedSeq2) const {
	//calculate the alignment score and record the traceback path
	const int numOfRows = seq1.length();
    // We operate on the full length of seq2 and build the full NW
    // matrix.
	const int numOfCols = seq2.length();
    // Setup the alignment and trace matrices
    const std::vector<int> dummyCol(numOfCols + 1);
    std::vector< std::vector<int> > alignMatrix(numOfRows + 1, dummyCol);
    // Setup direction matrices
    const std::vector<DIR> dummyDirEntries(numOfCols + 1);
    std::vector< std::vector<DIR> > trace(numOfRows + 1, dummyDirEntries);
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
        const int base1 = encodeBase(seq1[row - 1]); // Synonym
        // Compute the entries for current row and number of columns.
		for (int col = 1; (col <= numOfCols); col++) {
			// Initialize max to the first of the three terms (NORTH_WEST).
			const int base2 = encodeBase(seq2[col - 1]); // Synonym
			int maxScore    = scoreMatrix[base1][base2] +
                alignMatrix[row - 1][col - 1];
			DIR dirFlag  = NORTH_WEST;
			// See if the second term is larger (WEST).
			const int westScore = alignMatrix[row][col - 1] + gapPenalty;
			if (maxScore <= westScore) {
				maxScore = westScore;
				dirFlag  = WEST;
			}
			// See if the third term is the largest (NORTH)
			int northScore = alignMatrix[row - 1][col] + gapPenalty;
			if (maxScore <= northScore) {
				maxScore = northScore;
				dirFlag  = NORTH;
			}
            // Save the best score and direction of move in matrix for
            // use to reconstruct the aligned sequences below.
			alignMatrix[row][col] = maxScore;
			trace[row][col]       = dirFlag;
		}
	}
    std::for_each(alignMatrix.begin(), alignMatrix.end(), writeVector<int>);
    std::for_each(trace.begin(), trace.end(), writeVector<DIR>);
    // Now store the final alignment score in the parameter
	alignScore       = alignMatrix[numOfRows][numOfCols];
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
            ASSERT("Bad direction encountered in getNWAlignment()" == NULL);
            
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

int
AlignmentAlgorithm::getNWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   const int bandwidth,
                                   AlignmentResult& alignResult) const {
    return getNWAlignment(seq1, seq2, bandwidth, alignResult.alignScore,
                          alignResult.aScore, alignResult.alignedSeq1,
                          alignResult.alignedSeq2);
}

int
AlignmentAlgorithm::getNWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   const int bandwidth,
                                   int& alignScore, int& aScore,
                                   std::string& alignedSeq1,
                                   std::string& alignedSeq2) const {
	//calculate the alignment score and record the traceback path
	const int numOfRows = seq1.length() + 1;
    // If bandwidth is -1, then we don't use banded NW. We operate on
    // the full length of seq2 and build the full NW matrix.
	const int numOfCols = seq2.length() + 1;
    // Setup the alignment array that tracks only one row.
    std::vector<int> alignRow(numOfCols);
    // Setup direction matrices
    const std::vector<DIR> dummyDirEntries(bandwidth * 2);
    std::vector< std::vector<DIR> > trace(numOfRows, dummyDirEntries);
	trace[0][0] = BAD_DIR;
    alignRow[0] = 0;
	//initialize the first row and column of the score and trace
	//matrices as per needleman-wunch algorithm
	for (int col = 1; col < numOfCols; col++) {
		alignRow[col] = gapPenalty * col;
        trace[0][col] = WEST;
	}
	// Build the needleman-wunch (dynamic programming) matrix row-by-row
	for (int row = 1; row < numOfRows; row++) {
        int tmpMaxScore    = alignRow[0] + gapPenalty;
        const int base1    = encodeBase(seq1[row - 1]); // Synonym
        const int startCol = std::max(1, row - bandwidth);
        const int endCol   = std::min(row + bandwidth, numOfCols);
        const int colOff   = (row + bandwidth <= numOfCols) ? 1 :
            (bandwidth + bandwidth - (endCol - startCol));
        // Compute the entries for current row and diagonal.
		for (int col = startCol; (col < endCol); col++) {
			const int base2 = encodeBase(seq2[col - 1]); // Synonym
			// Initialize max to the first of the three terms (NORTH-WEST).
			int maxScore = alignRow[col - 1] + scoreMatrix[base1][base2];
			DIR dirFlag  = NORTH_WEST;
			// See if the second term is larger (WEST).
            if (std::abs(row - col + 1) <= bandwidth) {
                const int westScore = tmpMaxScore + gapPenalty;
                if (maxScore <= westScore) {
                    maxScore = westScore;
                    dirFlag  = WEST;
                }
            }
			// See if the third term is the largest (NORTH)
            if (abs(row - 1 - col) <= bandwidth) {
                const int northScore = alignRow[col] + gapPenalty;
                if (maxScore <= northScore) {
                    maxScore = northScore;
                    dirFlag  = NORTH;
                }
			}
            // Save the previous-maximum score for next iteration.
			alignRow[col - 1] = tmpMaxScore;
            // Setup the new maximum score for the next iteration
            tmpMaxScore = maxScore;
            // Save the best score and direction of move in matrix for
            // use to reconstruct the aligned sequences below.
			trace[row][col - startCol + colOff] = dirFlag;
            std::cout << dirFlag << " ";
		}
        std::cout << std::endl;
        alignRow[numOfCols - 1] = tmpMaxScore;
	}
    // Print the trace matrix for debugging purposes.
    std::for_each(trace.begin(), trace.end(), writeVector<DIR>);
    // Now store the final alignment score in the parameter
	alignScore       = alignRow[numOfCols - 1];
	// Trace back and get the alignment strings. Note that the
	// nucleotide sequence in tStr1 and tStr2 are reverse as we build
	// the alignments below.
    std::string tStr1;
    std::string tStr2;
    // Reserve sufficient memory in tStr1 and tStr2 to improve efficiency.
    const size_t maxLen = std::max(seq1.size(), seq2.size());
    tStr1.reserve(maxLen);
    tStr2.reserve(maxLen);

	int row    = numOfRows - 1;     // index into trace matrix & seq2
	int seqCol = numOfCols - 1;     // index into seq2
    int trCol  = bandwidth * 2 - 1; // index in trace matrix
    // Compute a-score as we go along constructing the aligned
    // sequences assuming seq1 is the reference sequence.
    aScore = 2 * seq1.size();
	while ((row > 0) || (trCol > 0)) {
		switch (trace[row][trCol]) {
		case NORTH:  // row-1, col, north
            row--;   // Must be done first as it affects code below
            trCol++; // Handle diagonal situation
			tStr1.append(1, seq1[row]);
			tStr2.append(1, '-');
            aScore += aScoreDelPen; // Deletion detected
			break;
		case WEST:    // row, seqCol-1, west
            seqCol--; // Must be done first as it affects code below
            trCol--;
			tStr1.append(1, '-');
			tStr2.append(1, seq2[seqCol]);
            aScore += aScoreInsPen; // Insertion detected
			break;
		case NORTH_WEST: // row-1, col-1, northwest
            row--;       // These two decrements must be done first as they
            seqCol--;    // affect seq1 and seq2 accesses right below.
			tStr1.append(1, seq1[row]);
			tStr2.append(1, seq2[seqCol]);
            // The nucleotides could be the same or different.
            aScore += (seq1[row] != seq2[seqCol] ? aScoreSubPen : 0);
			break;

        case BAD_DIR:
            ASSERT("Bad direction encountered in getNWAlignment()" == NULL);

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

int
AlignmentAlgorithm::getSWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   int& alignScore, int& aScore,
                                   std::string& alignedSeq1,
                                   std::string& alignedSeq2) const {
	//calculate the alignment score and record the traceback path
	const int numOfRows = seq1.length();
    // We operate on the full length of seq2 and build the full NW
    // matrix.
	const int numOfCols = seq2.length();
    // Setup the alignment and trace matrices
    const std::vector<int> dummyCol(numOfCols + 1);
    std::vector< std::vector<int> > alignMatrix(numOfRows + 1, dummyCol);
    // Setup direction matrices
    const std::vector<DIR> dummyDirEntries(numOfCols + 1);
    std::vector< std::vector<DIR> > trace(numOfRows + 1, dummyDirEntries);
	trace[0][0]       = BAD_DIR;
    alignMatrix[0][0] = 0;
	//initialize the first row and column of the score and trace
	//matrices as per smith-waterman algorithm
	for (int col = 1; col <= numOfCols; col++) {
		alignMatrix[0][col] = 0;
		trace[0][col]       = WEST;
	}
	for (int row = 1; row <= numOfRows; row++) {
		alignMatrix[row][0] = 0;
		trace[row][0]       = NORTH;
	}

	// Build the smith-waterman dynamic programming matrix row-by-row
    std::pair<int,int> maxScorePos;	//Location of Max Score in Matrix
    maxScorePos.first=0;
    maxScorePos.second=0;
	for (int row = 1; row <= numOfRows; row++) {
        const int base1 = encodeBase(seq1[row - 1]); // Synonym
        // Compute the entries for current row and number of columns.
		for (int col = 1; (col <= numOfCols); col++) {
			// Initialize max to the first of the three terms (NORTH_WEST).
			const int base2 = encodeBase(seq2[col - 1]); // Synonym
			int maxScore    = scoreMatrix[base1][base2] +
                alignMatrix[row - 1][col - 1];
			DIR dirFlag  = NORTH_WEST;
			// See if the second term is larger (WEST).
			const int westScore = alignMatrix[row][col - 1] + gapPenalty;
			if (maxScore <= westScore) {
				maxScore = westScore;
				dirFlag  = WEST;
			}
			// See if the third term is the largest (NORTH)
			int northScore = alignMatrix[row - 1][col] + gapPenalty;
			if (maxScore <= northScore) {
				maxScore = northScore;
				dirFlag  = NORTH;
			}
            // Ensure max score is never negative
            maxScore = std::min(0, maxScore);
			if (maxScore > alignMatrix[maxScorePos.first][maxScorePos.second]) {
				maxScorePos = std::make_pair(row, col);
			}
            // Save the best score and direction of move in matrix for
            // use to reconstruct the aligned sequences below.
			alignMatrix[row][col] = maxScore;
			trace[row][col]       = dirFlag;
		}
	}
    std::for_each(alignMatrix.begin(), alignMatrix.end(), writeVector<int>);
    std::for_each(trace.begin(), trace.end(), writeVector<DIR>);
    // Now store the final alignment score in the parameter
	alignScore       = alignMatrix[maxScorePos.first][maxScorePos.second];
	// Trace back and get the alignment strings. Note that the
	// nucleotide sequence in tStr1 and tStr2 are reverse as we build
	// the alignments below.
    std::string tStr1;
    std::string tStr2;
    // Reserve sufficient memory in tStr1 and tStr2 to improve efficiency.
    const size_t maxLen = std::max(seq1.size(), seq2.size());
    tStr1.reserve(maxLen);
    tStr2.reserve(maxLen);

	int row = maxScorePos.first;
	int col = maxScorePos.second;
    // Compute a-score as we go along constructing the aligned
    // sequences assuming seq1 is the reference sequence.
    aScore = 2 * seq1.size();
	while (((row > 0) || (col > 0) )&& alignMatrix[row][col]!=0) {
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
            ASSERT("Bad direction encountered in getNWAlignment()" == NULL);

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

int
AlignmentAlgorithm::getSWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   const int bandwidth,
                                   AlignmentResult& alignResult) const {
    return getSWAlignment(seq1, seq2, bandwidth, alignResult.alignScore,
                          alignResult.aScore, alignResult.alignedSeq1,
                          alignResult.alignedSeq2);
}

int
AlignmentAlgorithm::getSWAlignment(const std::string& seq1,
                                   const std::string& seq2,
                                   const int bandwidth,
                                   int& alignScore, int& aScore,
                                   std::string& alignedSeq1,
                                   std::string& alignedSeq2) const {
	//calculate the alignment score and record the traceback path
    const int numOfRows = seq1.length();
    // We operate on the full length of seq2 and build the full NW
    // matrix.
    const int numOfCols = bandwidth;
    // Setup the alignment and trace matrices
    const std::vector<int> dummyCol(numOfCols + 1);
    std::vector< std::vector<int> > alignMatrix(numOfRows + 1, dummyCol);
    // Setup direction matrices
    const std::vector<DIR> dummyDirEntries(numOfCols + 1);
    std::vector< std::vector<DIR> > trace(numOfRows + 1, dummyDirEntries);
    trace[0][0]       = BAD_DIR;
    alignMatrix[0][0] = 0;
    //initialize the first row and column of the score and trace
    //matrices as per smith-waterman algorithm
    for (int col = 1; col <= numOfCols; col++) {
        alignMatrix[0][col] = 0;
        trace[0][col]       = WEST;
    }
    for (int row = 1; row <= numOfRows; row++) {
        alignMatrix[row][0] = 0;
        trace[row][0]       = NORTH;
    }

    // Build the smith-waterman dynamic programming matrix row-by-row
    std::pair<int,int> maxScorePos;	//Location of Max Score in Matrix
    maxScorePos.first=0;
    maxScorePos.second=0;
    int colEnd=0;

    for (int row = 1; row <= numOfRows; row++) {
        const int base1 = encodeBase(seq1[row - 1]); // Synonym
        // Compute the entries for current row and number of columns equal to the bandwidth size
        for (int col = 1; (col <= numOfCols); col++) {
            // Initialize max to the first of the three terms (NORTH_WEST).
            const int base2 = encodeBase(seq2[col +colEnd-1]); // Synonym

            int maxScore    = scoreMatrix[base1][base2] +
                alignMatrix[row - 1][col];
            DIR dirFlag  = NORTH_WEST;
            // See if the second term is larger (WEST).
            const int westScore = alignMatrix[row][col - 1] + gapPenalty;
            if (maxScore <= westScore) {
                maxScore = westScore;
                dirFlag  = WEST;
            }
            // See if the third term is the largest (NORTH)
            int northScore = alignMatrix[row - 1][col+1] + gapPenalty;
            if (maxScore <= northScore) {
                maxScore = northScore;
                dirFlag  = NORTH;
            }
            if(maxScore<0){
                maxScore=0;
            }
            if(maxScore>alignMatrix[maxScorePos.first][maxScorePos.second]){
                maxScorePos.first=row;
                maxScorePos.second=col;
            }
            // Save the best score and direction of move in matrix for
            // use to reconstruct the aligned sequences below.
            alignMatrix[row][col] = maxScore;
            trace[row][col]       = dirFlag;
        }
        colEnd++;
    }
    std::for_each(alignMatrix.begin(), alignMatrix.end(), writeVector<int>);

    std::for_each(trace.begin(), trace.end(), writeVector<DIR>);
    // Now store the final alignment score in the parameter
    alignScore       = alignMatrix[maxScorePos.first][maxScorePos.second];
    // Trace back and get the alignment strings. Note that the
    // nucleotide sequence in tStr1 and tStr2 are reverse as we build
    // the alignments below.
    std::string tStr1;
    std::string tStr2;
    // Reserve sufficient memory in tStr1 and tStr2 to improve efficiency.
    const size_t maxLen = std::max(seq1.size(), seq2.size());
    tStr1.reserve(maxLen);
    tStr2.reserve(maxLen);

    int row = maxScorePos.first;
    int col = maxScorePos.second;
    // Compute a-score as we go along constructing the aligned
    // sequences assuming seq1 is the reference sequence.
    aScore = 2 * seq1.size();
    while (((row > 0) || (col > 0) )&& alignMatrix[row][col]!=0) {
        switch (trace[row][col]) {
        case NORTH: //row-1, col+1, north
            row--;
            col++;// Must be done first as it affects code below
            tStr1.append(1, seq1[row]);
            tStr2.append(1, '-');
            aScore += aScoreDelPen; // Deletion detected
            break;
        case WEST: //row, col-1, west
            col--; // Must be done first as it affects code below
            tStr1.append(1, '-');
            tStr2.append(1, seq2[col+row-1]);
            aScore += aScoreInsPen; // Insertion detected
            break;
        case NORTH_WEST: // row-1, col, northwest
            row--; // These two decrements msut be done first as they
            tStr1.append(1, seq1[row]);
            tStr2.append(1, seq2[col+row-1]);
            // The nucleotides could be the same or different.
            aScore += (seq1[row] != seq2[col+row-1] ? aScoreSubPen : 0);
            break;

        case BAD_DIR:
            ASSERT("Bad direction encountered in getNWAlignment()" == NULL);

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
