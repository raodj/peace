#include <math.h>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include "AlignmentAlgorithm.h"
using namespace std;

AlignmentAlgorithm::AlignmentAlgorithm(int match, int mismatch, int gap) {
	intMatrix score = createIntMatrix(6,6);
	score[0][0] =  match;
	score[0][1] =  mismatch;
	score[0][2] =  mismatch;
	score[0][3] =  mismatch;
	score[0][4] =  mismatch;
	score[0][5] =  mismatch;
	score[1][0] =  mismatch;
	score[1][1] =  match;
	score[1][2] =  mismatch;
	score[1][3] =  mismatch;
	score[1][4] =  mismatch;
	score[1][5] =  mismatch;
	score[2][0] =  mismatch;
	score[2][1] =  mismatch;
	score[2][2] =  match;
	score[2][3] =  mismatch;
	score[2][4] =  mismatch;
	score[2][5] =  mismatch;
	score[3][0] =  mismatch;
	score[3][1] =  mismatch;
	score[3][2] =  mismatch;
	score[3][3] =  match;
	score[3][4] =  mismatch;
	score[3][5] =  mismatch;
	score[4][0] =  mismatch;
	score[4][1] =  mismatch;
	score[4][2] =  mismatch;
	score[4][3] =  mismatch;
	score[4][4] =  match;
	score[4][5] =  mismatch;
	score[5][0] =  mismatch;
	score[5][1] =  mismatch;
	score[5][2] =  mismatch;
	score[5][3] =  mismatch;
	score[5][4] =  mismatch;
	score[5][5] =  match;

	this->setScoringMatrix(score);
	this->setGapPenalty(gap);

}

AlignmentAlgorithm::~AlignmentAlgorithm() {
	if (scoreMatrix != NULL)
		deleteIntMatrix(scoreMatrix,6);
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
	rows = s1.length();
	cols = s2.length();
	int encodedBases1[rows];
	int encodedBases2[cols];

	if (rows < cols) {
		// goes columnwise
		int array[rows];

		// initiate first column
		array[0] = 0;
		for (r = 1; r < rows; r++) {
			array[r] = array[r-1] + gapPenalty * r;
			encodedBases1[r-1] = encodeBase(s1[r-1]);
		}
		for (r=1; r < cols; r++)
			encodedBases2[r-1] = encodeBase(s2[r-1]);

		// calculate the similarity matrix (keep current column only)
		for (c = 1; c < cols; c++) {
			// initiate first row (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * c;
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
			array[c] = array[c-1] + gapPenalty * c;
			encodedBases2[c-1] = encodeBase(s2[c-1]);
		}
		for (r=1; r < rows; r++)
			encodedBases1[r-1] = encodeBase(s1[r-1]);

		// calculate the similarity matrix (keep current row only)
		for (r = 1; r < rows; r++)
		{
			// initiate first column (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * c;
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

int AlignmentAlgorithm::getBoundedNWScore(const string& s1, const string& s2) {
	int	r, c, rows, cols, tmp, ins, del, sub, score;
	int band = BAND_WIDTH_NW;
	rows = s1.length();
	cols = s2.length();
	int encodedBases1[rows];
	int encodedBases2[cols];

	if (rows < cols) {
		// goes columnwise
		int array[rows+1];

		// initiate first column
		array[0] = 0;
		for (r = 1; r <= rows; r++) {
			encodedBases1[r-1] = encodeBase(s1[r-1]);
			array[r] = array[r-1] + gapPenalty * r;
		}
		for (r=1; r <= cols; r++)
			encodedBases2[r-1] = encodeBase(s2[r-1]);

		// calculate the similarity matrix (keep current column only)
		for (c = 1; c <= cols; c++) {
			// initiate first row (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * c;
			int base2 = encodedBases2[c-1];
			int start = (c-band) > 1 ? (c-band) : 1;
			int end = cols > (c+band) ? (c+band) : cols;
			for (r = start; r < end; r++)
			{
				int base1 = encodedBases1[r-1];
				sub = array[r-1] + (this->scoreMatrix)[base1][base2];
				if (abs(c-r+1) <= band)
					ins = array[r] + gapPenalty;
				else
					ins = 0;
				if (abs(c-1-r) <= band)
					del = tmp + gapPenalty;
				else
					del = 0;

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
		for (c = 1; c <= cols; c++) {
			encodedBases2[c-1] = encodeBase(s2[c-1]);
			array[c] = array[c-1] + gapPenalty * c;
		}
		for (r=1; r <= rows; r++)
			encodedBases1[r-1] = encodeBase(s1[r-1]);

		// calculate the similarity matrix (keep current row only)
		for (r = 1; r <= rows; r++)
		{
			// initiate first column (tmp hold values
			// that will be later moved to the array)
			tmp = array[0] + gapPenalty * c;
			int base1 = encodedBases1[r-1];
			int start = (r-band) > 1 ? (r-band) : 1;
			int end = cols > (r+band) ? (r+band) : cols;

			for (c = start; c < end; c++)
			{
				int base2 = encodedBases2[c-1];
				sub = array[c-1] + (this->scoreMatrix)[base1][base2];
				if (abs(c-1-r) <= band)
					ins = tmp + gapPenalty;
				else
					ins = 0;
				if (abs(r-c+1) <= band)
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
	}
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
	//intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	int trace[numOfRows + 1][numOfCols + 1];
	trace[0][0] = 0;

	//intMatrix alignMatrix = createIntMatrix(numOfRows + 1, numOfCols + 1);
	int alignMatrix[numOfRows + 1][numOfCols + 1];
	//initialize the matrix
	alignMatrix[0][0] = 0;
	//int* encodedBases2 = new int[numOfCols];
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

	//deleteIntMatrix(trace, numOfRows+1);
	//deleteIntMatrix(alignMatrix, numOfRows+1);
	//delete [] encodedBases2;

	return result;
}

AlignResult AlignmentAlgorithm::getBoundedSWAlignment(const std::string& horizontal,
		const std::string& vertical) {
	AlignResult result;

	//calculate the alignment score and record the traceback path
	int numOfRows = vertical.length();
	int numOfCols = horizontal.length();
	//intMatrix trace = createIntMatrix(numOfRows + 1, numOfCols + 1);
	int trace[numOfRows + 1][numOfCols + 1];
	trace[0][0] = 0;

	int alignMatrix[numOfRows + 1][numOfCols + 1];
	//initialize the matrix
	alignMatrix[0][0] = 0;
	int encodedBaseH[numOfCols];
	for (int i = 1; i <= numOfCols; i++) {
		alignMatrix[0][i] = 0;
		trace[0][i] = 0; //start point
		encodedBaseH[i-1] = encodeBase(horizontal[i-1]);
	}

	int maxScore = 0;
	int maxRow = 0;
	int maxCol = 0;

	//build the matrix row by row until row 50 and then look for the starting point of bounded algorithm in the first 50 rows.
	int tmpRowEnd = numOfRows>50? 50 : numOfRows;
	for (int i = 1; i <= tmpRowEnd; i++) {
		int base1 = encodeBase(vertical[i - 1]);
		alignMatrix[i][0] = 0;
		trace[i][0] = 0; //start point

		for (int j = 1; j <= numOfCols; j++) {
			int flag = 1;
			// Initialize max to the first of the three terms (NORTH).
			int base2 = encodedBaseH[j - 1];
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

	//check trace to decide where to start bounded algorithm
	int boundRow = 0;
	int boundCol = 0;
	int bounded = 0;
	int band = BAND_WIDTH_SW;
	for (int i=0; i<=tmpRowEnd-band; i++) {
		for (int j=0; j<=numOfCols-band; j++) {
			if (alignMatrix[i][j] == 0) {
				int northwestNum = 0;
				for (int k=1; k<=band; k++)
					northwestNum = trace[i+k][j+k]==3 ? (northwestNum+1): northwestNum;

				if (northwestNum >= band-2) {
					bounded = 1;
					boundRow = i;
					boundCol = j;
					break;
				}
			}
		}
		if (bounded == 1) break;
	}

	//build the matrix row by row
	int tmpStart = boundCol+1;
	for (int i = boundRow+1; i <= numOfRows; i++) {
		int base1 = encodeBase(vertical[i - 1]);

		int tmpI = i + boundCol;
		int start = (tmpI-band) > tmpStart ? (tmpI-band) : tmpStart;
		int end = numOfCols > (tmpI+band) ? (tmpI+band) : numOfCols;
		for (int j = start; j <= end; j++) {
			int base2 = encodedBaseH[j - 1];
			int flag = 3;
			// Initialize max to the first of the three terms (NORTHWEST).
			int max = alignMatrix[i - 1][j - 1]
		                                   + (this->scoreMatrix)[base1][base2];

			// See if the second term is larger (WEST).
			if (abs(tmpI-j+1) <= band) {
				int west = alignMatrix[i][j - 1] + this->gapPenalty;
				if (max <= west) {
					max = west;
					flag = 2;
				}
			}

			// See if the third term is the largest (NORTH)
			if (abs(tmpI-1-j) <= band) {
				int north = alignMatrix[i - 1][j] + this->gapPenalty;
				if (max <= north) {
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
			tStr1.append(1, vertical[row-1]);
			tStr2.append(1, '-');
			row = row - 1;
			break;
		case 2: //i, j-1, west
			tStr1.append(1, '-');
			tStr2.append(1, horizontal[col-1]);
			col = col - 1;
			break;
		case 3: //i-1, j-1, northwest
			tStr1.append(1, vertical[row-1]);
			tStr2.append(1, horizontal[col-1]);
			row = row - 1;
			col = col - 1;
			break;
		}
		flag = trace[row][col];
	}

	//set str1 and str2 in result, they are reverse of tStr1 and tStr2
	reverse(tStr1.begin(), tStr1.end());
	reverse(tStr2.begin(), tStr2.end());
	result.str1 = tStr2; //result.str1 for the first parameter "horizontal"
	result.str2 = tStr1; //result.str2 for the second parameter "vertical"
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

