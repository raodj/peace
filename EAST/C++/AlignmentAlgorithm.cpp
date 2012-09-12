#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include "AlignmentAlgorithm.h"
using namespace std;

AlignmentAlgorithm::AlignmentAlgorithm(int match, int mismatch, int gap) {
	this->setScoringMatrix(match, mismatch);
	this->setGapPenalty(gap);

}

AlignmentAlgorithm::~AlignmentAlgorithm() {
	if (scoreMatrix != NULL)
		deleteIntMatrix(scoreMatrix,6);
}

void AlignmentAlgorithm::setScoringMatrix(int match, int mismatch) {
	scoreMatrix = createIntMatrix(17, 17);
	int m = match;
	int n = mismatch;
					   //A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X
	int score[17][17] = {{m, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n},  //A
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

	for (int i=0; i<17; i++) {
		for (int j=0; j<17; j++) {
			scoreMatrix[i][j] = score[i][j];
		}
	}
}

int AlignmentAlgorithm::getNWScore(const string& s1, const string& s2) {
	/*
	int numOfRows = s1.length();
	int numOfCols = s2.length();

	intMatrix alignMatrix = createIntMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	alignMatrix[0][0] = 0;
	int* encodedBases2 = new int[numOfCols];
	for (int i = 1; i <= numOfCols; i++) {
	  alignMatrix[0][i] = gapPenalty * i;
	  encodedBases2[i-1] = encodeBase(s2[i-1]);
	}

	//build the matrix row by row
	int base1, base2, north, west, northwest;
	for (int i = 1; i <= numOfRows; i++) {
	  base1 = encodeBase(s1[i-1]);
	  // initiate first column
	  alignMatrix[i][0] = gapPenalty * i;

	  for (int j = 1; j <= numOfCols; j++) {
	    // Initialize max to the first of the three terms (NORTH).
	    base2 = encodedBases2[j - 1];
	    north = alignMatrix[i - 1][j] + gapPenalty;

	    // See if the second term is larger (WEST).
	    west = alignMatrix[i][j - 1] + gapPenalty;

	    // See if the third term is the largest (NORTHWEST)
	    northwest = alignMatrix[i - 1][j - 1]
	      + (this->scoreMatrix)[base1][base2];

	    if (northwest >= north)
	      alignMatrix[i][j] = northwest >= west ? northwest : west;
	    else
	      alignMatrix[i][j] = north >= west ? north : west;
	  }
	}
	int score = alignMatrix[numOfRows][numOfCols];

	deleteIntMatrix(alignMatrix, numOfRows+1);
	delete [] encodedBases2;
	return score;
	 */
	int		r, c, rows, cols, tmp, ins, del, sub, score;
	rows = s1.length() + 1;
	cols = s2.length() + 1;
	int encodedBases1[rows];
	int encodedBases2[cols];

	if (rows < cols) {
		// goes columnwise
		int array[rows];

		// initiate first column
		array[0] = 0;
		for (r = 1; r < rows; r++) {
			array[r] = array[0] + gapPenalty * r;
			encodedBases1[r-1] = encodeBase(s1[r-1]);
		}
		for (r=1; r < cols; r++)
			encodedBases2[r-1] = encodeBase(s2[r-1]);

		// calculate the similarity matrix (keep current column only)
		for (c = 1; c < cols; c++) {
			// initiate first row (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty;
			int base2 = encodedBases2[c-1];
			for (r = 1; r < rows; r++)
			{
				int base1 = encodedBases1[r-1];
				ins = array[r] + gapPenalty;
				sub = array[r-1] + (this->scoreMatrix)[base1][base2];
				del = tmp + gapPenalty;

				// move the temp value to the array
				array[r-1] = tmp;

				// choose the greatest
				if (sub >= ins)
					tmp = sub >= del ? sub : del;
				else
					tmp = ins >= del ? ins : del;
			}

			// move the temp value to the array
			array[rows - 1] = tmp;
		}
		score = array[rows - 1];
	} else {
		// goes rowwise
		int array[cols];

		// initiate first row
		array[0] = 0;
		for (c = 1; c < cols; c++) {
			array[c] = array[0] + gapPenalty * c;
			encodedBases2[c-1] = encodeBase(s2[c-1]);
		}
		for (r=1; r < rows; r++)
			encodedBases1[r-1] = encodeBase(s1[r-1]);

		// calculate the similarity matrix (keep current row only)
		for (r = 1; r < rows; r++)
		{
			// initiate first column (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty;
			int base1 = encodedBases1[r-1];

			for (c = 1; c < cols; c++)
			{
				int base2 = encodedBases2[c-1];
				ins = tmp + gapPenalty;
				sub = array[c-1] + (this->scoreMatrix)[base1][base2];
				del = array[c] + gapPenalty;

				// move the temp value to the array
				array[c-1] = tmp;

				// choose the greatest
				if (sub >= ins)
					tmp = sub >= del ? sub : del;
				else
					tmp = ins >= del ? ins : del;
			}

			// move the temp value to the array
			array[cols - 1] = tmp;
		}
		score = array[cols - 1];
	}
	return score;
}

int AlignmentAlgorithm::getNWScoreWithQual(const string& s1, const string& s2, const std::vector<int>& qualScore1, const std::vector<int>& qualScore2) {
	int		r, c, rows, cols, tmp, ins, del, sub, score;
	rows = s1.length() + 1;
	cols = s2.length() + 1;
	int encodedBases1[rows];
	int encodedBases2[cols];

	if (rows < cols) {
		// goes columnwise
		int array[rows];

		// initiate first column
		array[0] = 0;
		for (r = 1; r < rows; r++) {
			array[r] = array[r-1] + gapPenalty * qualScore1[r-1];
			encodedBases1[r-1] = encodeBase(s1[r-1]);
		}
		for (r=1; r < cols; r++)
			encodedBases2[r-1] = encodeBase(s2[r-1]);

		// calculate the similarity matrix (keep current column only)
		for (c = 1; c < cols; c++) {
			// initiate first row (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * qualScore2[c-1];
			int base2 = encodedBases2[c-1];
			for (r = 1; r < rows; r++)
			{
				int base1 = encodedBases1[r-1];
				int score1 = 0;
				int score2 = 0;
				if (r < rows-1) {
					score1 = min(qualScore1[r-1], qualScore1[r]);
				} else {
					score1 = qualScore1[r-1];
				}
				if (c < cols-1) {
					score2 = min(qualScore2[c-1], qualScore2[c]);
				} else {
					score2 = qualScore2[c-1];
				}
				ins = array[r] + gapPenalty * min(qualScore2[c-1], score1);
				sub = (array[r-1] + (this->scoreMatrix)[base1][base2]) * min(qualScore1[r-1], qualScore2[c-1]);
				del = tmp + gapPenalty * min(qualScore1[r-1], score2);

				// move the temp value to the array
				array[r-1] = tmp;

				// choose the greatest
				if (sub >= ins)
					tmp = sub >= del ? sub : del;
				else
					tmp = ins >= del ? ins : del;
			}

			// move the temp value to the array
			array[rows - 1] = tmp;
		}
		score = array[rows - 1];
	} else {
		// goes rowwise
		int array[cols];

		// initiate first row
		array[0] = 0;
		for (c = 1; c < cols; c++) {
			array[c] = array[c-1] + gapPenalty * qualScore2[c-1];
			encodedBases2[c-1] = encodeBase(s2[c-1]);
		}
		for (r=1; r < rows; r++)
			encodedBases1[r-1] = encodeBase(s1[r-1]);

		// calculate the similarity matrix (keep current row only)
		for (r = 1; r < rows; r++)
		{
			// initiate first column (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * qualScore1[r-1];
			int base1 = encodedBases1[r-1];

			for (c = 1; c < cols; c++)
			{
				int base2 = encodedBases2[c-1];
				int score1 = 0;
				int score2 = 0;
				if (r < rows-1) {
					score1 = min(qualScore1[r-1], qualScore1[r]);
				} else {
					score1 = qualScore1[r-1];
				}
				if (c < cols-1) {
					score2 = min(qualScore2[c-1], qualScore2[c]);
				} else {
					score2 = qualScore2[c-1];
				}

				ins = tmp + gapPenalty * min(qualScore2[c-1], score1);
				sub = array[c-1] + (this->scoreMatrix)[base1][base2] * min(qualScore1[r-1], qualScore2[c-1]);
				del = array[c] + gapPenalty * min(qualScore1[r-1], score2);

				// move the temp value to the array
				array[c-1] = tmp;

				// choose the greatest
				if (sub >= ins)
					tmp = sub >= del ? sub : del;
				else
					tmp = ins >= del ? ins : del;
			}

			// move the temp value to the array
			array[cols - 1] = tmp;
		}
		score = array[cols - 1];
	}
	return score;
}

int AlignmentAlgorithm::getBoundedNWScore(const string& s1, const string& s2) {
	int	r, c, rows, cols, tmp, ins, del, sub, score;
	int band = BAND_WIDTH_NW;
	string ss1;
	string ss2;
	if (s1.length() < s2.length()) {
		ss1 = s2;
		ss2 = s1;
	} else {
		ss1 = s1;
		ss2 = s2;
	}
	rows = ss1.length() + 1;
	cols = ss2.length() + 1;
	int encodedBases1[rows-1];
	int encodedBases2[cols-1];

	// goes rowwise
	int array[cols];

	// initiate first row
	array[0] = 0;
	for (c = 1; c < cols; c++) {
		encodedBases2[c-1] = encodeBase(ss2[c-1]);
		array[c] = array[0] + gapPenalty * c;
	}
	for (r=1; r < rows; r++)
		encodedBases1[r-1] = encodeBase(ss1[r-1]);

	// calculate the similarity matrix (keep current row only)
	for (r = 1; r < rows; r++)
	{
		tmp = array[0] + gapPenalty;
		int base1 = encodedBases1[r-1];
		int start = (r-band) > 1 ? (r-band) : 1;
		int end = cols > (r+band) ? (r+band) : cols;

		for (c = start; c < end; c++)
		{
			int base2 = encodedBases2[c-1];
			sub = array[c-1] + (this->scoreMatrix)[base1][base2];
			if (abs(r-c+1) <= band) //west
				ins = tmp + gapPenalty;
			else
				ins = 0;
			if (abs(c-1-r) <= band) //north
				del = array[c] + gapPenalty;
			else
				del = 0;

			// move the temp value to the array
			array[c-1] = tmp;

			// choose the greatest
			if (sub >= ins)
				tmp = sub >= del ? sub : del;
			else
				tmp = ins >= del ? ins : del;
		}

		// move the temp value to the array
		array[cols - 1] = tmp;
	}
	score = array[cols - 1];

	return score;
}

int AlignmentAlgorithm::getBoundedNWScoreWithQual(const string& s1, const string& s2, const std::vector<int>& qual1, const std::vector<int>& qual2) {
	int	r, c, rows, cols, tmp, ins, del, sub, score;
	int band = BAND_WIDTH_NW;
	string ss1;
	string ss2;
	vector<int> qualScore1;
	vector<int> qualScore2;
	if (s1.length() < s2.length()) {
		ss1 = s2;
		ss2 = s1;
		qualScore1 = qual2;
		qualScore2 = qual1;
	} else {
		ss1 = s1;
		ss2 = s2;
		qualScore1 = qual1;
		qualScore2 = qual2;
	}
	rows = ss1.length() + 1;
	cols = ss2.length() + 1;
	int encodedBases1[rows-1];
	int encodedBases2[cols-1];

	// goes rowwise
	int array[cols];

	// initiate first row
	array[0] = 0;
	for (c = 1; c < cols; c++) {
		encodedBases2[c-1] = encodeBase(ss2[c-1]);
		array[c] = array[c-1] + gapPenalty * qualScore2[c-1];
	}
	for (r=1; r < rows; r++)
		encodedBases1[r-1] = encodeBase(ss1[r-1]);

	// calculate the similarity matrix (keep current row only)
	for (r = 1; r < rows; r++)
	{
		tmp = array[0] + gapPenalty * qualScore1[r-1];
		int base1 = encodedBases1[r-1];
		int start = (r-band) > 1 ? (r-band) : 1;
		int end = cols > (r+band) ? (r+band) : cols;

		for (c = start; c < end; c++)
		{
			int base2 = encodedBases2[c-1];
			sub = array[c-1] + (this->scoreMatrix)[base1][base2] * min(qualScore1[r-1], qualScore2[c-1]);

			int score1 = 0;
			int score2 = 0;
			if (r < rows-1) {
				score1 = min(qualScore1[r-1], qualScore1[r]);
			} else {
				score1 = qualScore1[r-1];
			}
			if (c < cols-1) {
				score2 = min(qualScore2[c-1], qualScore2[c]);
			} else {
				score2 = qualScore2[c-1];
			}

			if (abs(r-c+1) <= band) //west
				ins = tmp + gapPenalty * min(qualScore2[c-1], score1);
			else
				ins = 0;
			if (abs(c-1-r) <= band) //north
				del = array[c] + gapPenalty * min(qualScore1[r-1], score2);
			else
				del = 0;

			// move the temp value to the array
			array[c-1] = tmp;

			// choose the greatest
			if (sub >= ins)
				tmp = sub >= del ? sub : del;
			else
				tmp = ins >= del ? ins : del;
		}

		// move the temp value to the array
		array[cols - 1] = tmp;
	}
	score = array[cols - 1];

	return score;
}

AlignResult AlignmentAlgorithm::getNWAlignment(const std::string& s1,
		const std::string& s2) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();
	intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	trace[0][0] = 0;

	intMatrix alignMatrix = createIntMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix[0][i] = (this->gapPenalty) * i;
		trace[0][i] = 2; //west
	}
	for (int i = 1; i <= numOfRows; i++) {
		alignMatrix[i][0] = (this->gapPenalty) * i;
		trace[i][0] = 1; //north
	}

	vector<int> encodedBases1;
	vector<int> encodedBases2;
	for (int i=0; i<=numOfRows-1; i++) {
		encodedBases1.push_back(encodeBase(s1[i]));
	}
	for (int i=0; i<=numOfCols-1; i++) {
		encodedBases2.push_back(encodeBase(s2[i]));
	}
	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodedBases1[i - 1];
			int base2 = encodedBases2[j - 1];
			int max = alignMatrix[i - 1][j] + this->gapPenalty;

			// See if the second term is larger (WEST).
			int west = alignMatrix[i][j - 1] + this->gapPenalty;
			if (max <= west) {
				max = west;
				flag = 2;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix[i - 1][j - 1]
			                                   + (this->scoreMatrix)[base1][base2];
			if (max <= northwest) {
				max = northwest;
				flag = 3;
			}

			alignMatrix[i][j] = max;
			trace[i][j] = flag;
		}
	}

	result.score = alignMatrix[numOfRows][numOfCols];

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = numOfRows;
	int col = numOfCols;
	while ((row != 0) || (col != 0)) {
		int flag = trace[row][col];
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	reverse(tStr1.begin(), tStr1.end());
	reverse(tStr2.begin(), tStr2.end());
	result.str1 = tStr1;
	result.str2 = tStr2;

	deleteIntMatrix(trace, numOfRows+1);
	deleteIntMatrix(alignMatrix, numOfRows+1);

	return result;
}

int AlignmentAlgorithm::getSWScore(const string& s1, const string& s2) {
	int numOfRows = s1.length();
	int numOfCols = s2.length();

	intMatrix alignMatrix = createIntMatrix(numOfRows + 1, numOfCols + 1);

	//initialize the matrix
	for (int i = 0; i <= numOfCols; i++) {
		alignMatrix[0][i] = 0;
	}
	for (int i = 0; i <= numOfRows; i++) {
		alignMatrix[i][0] = 0;
	}

	int maxScore = 0;
	vector<int> encodedBases1;
	vector<int> encodedBases2;
	for (int i=0; i<=numOfRows-1; i++) {
		encodedBases1.push_back(encodeBase(s1[i]));
	}
	for (int i=0; i<=numOfCols-1; i++) {
		encodedBases2.push_back(encodeBase(s2[i]));
	}
	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		for (int j = 1; j <= numOfCols; j++) {
			// Initialize max to the first of the three terms (NORTH).
			int base1 = encodedBases1[i - 1];
			int base2 = encodedBases2[j - 1];
			int innerMax = max(0, alignMatrix[i - 1][j] + this->gapPenalty);

			// See if the second term is larger (WEST).
			int west = alignMatrix[i][j - 1] + this->gapPenalty;
			if (innerMax <= west) {
				innerMax = west;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix[i - 1][j - 1]
			                                   + (this->scoreMatrix)[base1][base2];
			if (innerMax <= northwest) {
				innerMax = northwest;
			}

			alignMatrix[i][j] = innerMax;

			maxScore = max(maxScore, innerMax);
		}
	}

	deleteIntMatrix(alignMatrix, numOfRows+1);
	return maxScore;
}

AlignResult AlignmentAlgorithm::getSWAlignment(const std::string& s1,
		const std::string& s2) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();

	intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	intMatrix alignMatrix = createIntMatrix(numOfRows + 1, numOfCols + 1);

	trace[0][0] = 0;
	alignMatrix[0][0] = 0;
    
	//int* encodedBases2 = new int[numOfCols];
    std::vector<int> encodedBases2(numOfCols);
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix[0][i] = 0;
		trace[0][i] = 0; //start point
		encodedBases2[i-1] = encodeBase(s2[i-1]);
	}

	int maxScore = 0;
	int maxRow = 0;
	int maxCol = 0;

	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		int base1 = encodeBase(s1[i - 1]);
		alignMatrix[i][0] = 0;
		trace[i][0] = 0; //start point

		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			// Initialize max to the first of the three terms (NORTH).
			int base2 = encodedBases2[j - 1];
			int max = alignMatrix[i - 1][j] + this->gapPenalty;

			// See if the second term is larger (WEST).
			int west = alignMatrix[i][j - 1] + this->gapPenalty;
			if (max <= west) {
				max = west;
				flag = 2;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix[i - 1][j - 1]
			                                   + (this->scoreMatrix)[base1][base2];
			if (max <= northwest) {
				max = northwest;
				flag = 3;
			}

			if (max <= 0) {
				alignMatrix[i][j] = 0;
				trace[i][j] = 0; //start point
			} else {
				alignMatrix[i][j] = max;
				trace[i][j] = flag;
			}
			if (max > maxScore) {
				maxScore = max;
				maxRow = i;
				maxCol = j;
			}
		}
	}

	result.score = alignMatrix[maxRow][maxCol];

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = maxRow;
	int col = maxCol;
	int flag = trace[row][col];
	while (flag != 0) {
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
		flag = trace[row][col];
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	reverse(tStr1.begin(), tStr1.end());
	reverse(tStr2.begin(), tStr2.end());
	result.str1 = tStr1;
	result.str2 = tStr2;

	deleteIntMatrix(trace, numOfRows+1);
	deleteIntMatrix(alignMatrix, numOfRows+1);

	//delete [] encodedBases2;

	return result;
}

AlignResult AlignmentAlgorithm::getSWAlignmentWithQual(const std::string& s1,
		const std::string& s2, const std::vector<int>& qualScore1, const std::vector<int>& qualScore2) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();
	//intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	int trace[numOfRows + 1][numOfCols + 1];
	trace[0][0] = 0;

	int alignMatrix[numOfRows + 1][numOfCols + 1];
	//initialize the matrix
	alignMatrix[0][0] = 0;
	int encodedBases2[numOfCols];
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix[0][i] = 0;
		trace[0][i] = 0; //start point
		encodedBases2[i-1] = encodeBase(s2[i-1]);
	}

	int maxScore = 0;
	int maxRow = 0;
	int maxCol = 0;

	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		int base1 = encodeBase(s1[i - 1]);
		alignMatrix[i][0] = 0;
		trace[i][0] = 0; //start point

		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			int score1 = 0;
			int score2 = 0;
			if (i < numOfRows) {
				score1 = min(qualScore1[i-1], qualScore1[i]);
			} else {
				score1 = qualScore1[i-1];
			}
			if (j < numOfCols) {
				score2 = min(qualScore2[j-1], qualScore2[j]);
			} else {
				score2 = qualScore2[j-1];
			}

			// Initialize max to the first of the three terms (NORTH).
			int base2 = encodedBases2[j - 1];
			int max = alignMatrix[i - 1][j] + this->gapPenalty * min(qualScore1[i-1], score2);

			// See if the second term is larger (WEST).
			int west = alignMatrix[i][j - 1] + this->gapPenalty * min(qualScore2[j-1], score1);
			if (max <= west) {
				max = west;
				flag = 2;
			}

			// See if the third term is the largest (NORTHWEST)
			int northwest = alignMatrix[i - 1][j - 1]
			                                   + (this->scoreMatrix)[base1][base2] * min(qualScore1[i-1], qualScore2[j-1]);
			if (max <= northwest) {
				max = northwest;
				flag = 3;
			}

			if (max <= 0) {
				alignMatrix[i][j] = 0;
				trace[i][j] = 0; //start point
			} else {
				alignMatrix[i][j] = max;
				trace[i][j] = flag;
			}
			if (max > maxScore) {
				maxScore = max;
				maxRow = i;
				maxCol = j;
			}
		}
	}

	result.score = alignMatrix[maxRow][maxCol];

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = maxRow;
	int col = maxCol;
	int flag = trace[row][col];
	while (flag != 0) {
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
		flag = trace[row][col];
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	reverse(tStr1.begin(), tStr1.end());
	reverse(tStr2.begin(), tStr2.end());
	result.str1 = tStr1;
	result.str2 = tStr2;

	return result;
}

AlignResult AlignmentAlgorithm::getBoundedSWAlignment(const std::string& s1,
		const std::string& s2) {
	AlignResult result;
	int band = BAND_WIDTH_SW;

	//calculate the alignment score and record the traceback path
	int numOfRows = s1.length();
	int numOfCols = s2.length();
	//intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	int trace[numOfRows + 1][numOfCols + 1];
	trace[0][0] = 0;

	int alignMatrix[numOfRows + 1][numOfCols + 1];
	//initialize the matrix
	alignMatrix[0][0] = 0;
	int encodedBase2[numOfCols];
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix[0][i] = 0;
		trace[0][i] = 0; //start point
		encodedBase2[i-1] = encodeBase(s2[i-1]);
	}

	for (int i=0; i<=numOfRows; i++)
		for (int j=0; j<=numOfCols; j++)
			alignMatrix[i][j] = -100;
	alignMatrix[0][0] = 0;

	int maxScore = 0;
	int maxRow = 0;
	int maxCol = 0;


	//build the matrix row by row
	for (int i = 1; i <= numOfRows; i++) {
		int base1 = encodeBase(s1[i - 1]);
		alignMatrix[i][0] = 0;
		trace[i][0] = 0;

		int start = (i-band) > 1 ? (i-band) : 1;
		int end = numOfCols > (i+band) ? (i+band) : numOfCols;
		for (int j = start; j <= end; j++) {
			int base2 = encodedBase2[j - 1];
			int flag = 3;
			// Initialize max to the first of the three terms (NORTHWEST).
			int max = alignMatrix[i - 1][j - 1]
		                                   + (this->scoreMatrix)[base1][base2];

			// See if the second term is larger (WEST).
			if (abs(i-j+1) <= band) {
				int west = alignMatrix[i][j - 1] + this->gapPenalty;
				if (max < west) {
					max = west;
					flag = 2;
				}
			}

			// See if the third term is the largest (NORTH)
			if (abs(i-1-j) <= band) {
				int north = alignMatrix[i - 1][j] + this->gapPenalty;
				if (max < north) {
					max = north;
					flag = 1;
				}
			}

			if (max <= 0) {
				alignMatrix[i][j] = 0;
				trace[i][j] = 0; //start point
			} else {
				alignMatrix[i][j] = max;
				trace[i][j] = flag;
			}
			if (max > maxScore) {
				maxScore = max;
				maxRow = i;
				maxCol = j;
			}
		}
	}

	result.score = alignMatrix[maxRow][maxCol];

	//trace back and get the alignment strings
	string tStr1;
	string tStr2;
	int row = maxRow;
	int col = maxCol;
	int flag = trace[row][col];
	while (flag != 0) {
		switch (flag) {
		case 1: //i-1, j, north
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, s2[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, s1[row-1]);
			tStr2.append(1, s2[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
		flag = trace[row][col];
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	reverse(tStr1.begin(), tStr1.end());
	reverse(tStr2.begin(), tStr2.end());
	result.str1 = tStr1;
	result.str2 = tStr2;
	return result;
}

intMatrix createIntMatrix(int n, int m) {
	int** p = new int*[n];
	int** p2;
	int i;
	for (i=0, p2=p; i < n; i++, p2++) {
		*p2  = new int[m];
	}
	return p;
}
void deleteIntMatrix(int** p, int n) {
	int** p2;
	int i;
	for (i=0, p2=p; i < n; i++, p2++)
		delete [] *p2;
	delete [] p;
}


