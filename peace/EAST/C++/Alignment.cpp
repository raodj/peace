#include "Alignment.h"
#include <time.h>

using namespace std;

Alignment::Alignment() {
	alignmentThreshold = ALIGNMENT_THRESHOLD;
	alignAlgo = new AlignmentAlgorithm(1,-1,-2);
	numCall = 0;
	usedTime = 0;
	numCall2 = 0;
	usedTime2 = 0;
}

Alignment::~Alignment() {
  delete alignAlgo;
}

/*
 * Calculate the distance of two strings.
 * dis = (1 - similarityScore/lengthOfLongerString)*a, actually, in our case, s1 has the same length as s2.
 * Now, we set a=100. So the return value would be [0, 100]
 * @param s1, s2
 * @return int distance.
 */
int Alignment::getDistance(const string& s1, const string& s2, vector<int> qualScore1, vector<int> qualScore2) {
	this->numCall2++;
	time_t start,end;
	start = time(NULL);
	int score;
	if (USE_QUALITY_FILE == 1) { //use quality file
		score = getSimlarityScoreWithQual(s1, s2, qualScore1, qualScore2);
	} else {
		score = getSimlarityScore(s1, s2);
	}

	end = time(NULL);
	usedTime2 += end - start;

	int retVal = INT_MAX;
	if (score != 0) {
		int length = s1.length();
		if (USE_QUALITY_FILE == 1) { //use quality file
			double totalQual = 0;
			for (int i=0; i<length; i++) {
				totalQual += qualScore1[i] + qualScore2[i];
			}
			retVal = (int) ((1 - (double) score / (totalQual/2)) * 100);
			if (retVal < 0) {
				retVal = 0;
			}
		} else {
			retVal = (int) ((1 - (double) score / length) * 100);
		}
	}

	if (retVal > alignmentThreshold) {
		retVal = INT_MAX;
	}
	return retVal;
}

/*
 * Use Needleman-Wunsch algorithm to calculate similarity score of two string.
 * @param s1, s2
 * @return int similarity score(>=0), if the value is less than 0, it's set to be zero.
 */
int Alignment::getSimlarityScore(const string& s1, const string& s2) {
	if (USE_BOUNDED_NW == 0) { //use ordinary version
		return alignAlgo->getNWScore(s1, s2);
	} else { //use bounded version
		return alignAlgo->getBoundedNWScore(s1, s2);
	}
}

int Alignment::getSimlarityScoreWithQual(const string& s1, const string& s2, const vector<int>& qualScore1, const vector<int>& qualScore2) {
	if (USE_BOUNDED_NW == 0) { //use ordinary version
		return alignAlgo->getNWScoreWithQual(s1, s2, qualScore1, qualScore2);
	} else { //use bounded version
		return alignAlgo->getBoundedNWScoreWithQual(s1, s2, qualScore1, qualScore2);
	}
}

/*
 * Use Smith-Waterman algorithm to get local alignment.
 * @param s1, s2
 * @return string[], [0] and [1] are the two aligned sequences, [2] is the pairwise alignment.
 */
AlignResult Alignment::getLocalAlignment(const string& s1, const string& s2, vector<int> qualScore1, vector<int> qualScore2) {
	//return alignAlgo->getSWAlignment(s1, s2);
	this->numCall++;
	time_t start,end;
	start = time(NULL);
	AlignResult ret;
	if (USE_QUALITY_FILE == 1) { //use quality file
		ret = alignAlgo->getSWAlignmentWithQual(s1, s2, qualScore1, qualScore2);
	} else {
		ret = alignAlgo->getSWAlignment(s1, s2);
	}
	end = time(NULL);
	usedTime += end - start;
	return ret;
}

/*
 * Use Bounded Smith-Waterman algorithm to get local alignment. This method is used to speed up.
 * @param s1, s2
 * @return string[], [0] and [1] are the two aligned sequences, [2] is the pairwise alignment.
 */
AlignResult Alignment::getBoundedLocalAlignment(const string& s1, const string& s2) {
	//return alignAlgo->getBoundedSWAlignment(s1, s2);
	this->numCall++;
	time_t start,end;
	start = time(NULL);
	AlignResult ret = alignAlgo->getBoundedSWAlignment(s1, s2);
	end = time(NULL);
	usedTime += end - start;
	return ret;
}

void Alignment::setScoringSystem(int match, int mismatch, int gap) {
	alignAlgo->setScoringMatrix(match, mismatch);
	alignAlgo->setGapPenalty(gap);
}
