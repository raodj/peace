#include "OvlDistance.h"
#include <cstdlib>

using namespace std;

/**
 * Get the overlap distance of the two strings.
 * This function
 * 	1) tries to find the position with the minimal d2 distance;
 * 	2) uses alignment to get the similarity value of two substrings. And it sets
 * 		the similarity value to be the overlap distance of s1 and s2.
 *
 * The function finds the position in s2 which has the smallest d2 distance between s1's first window and
 * s2 or between s1's last window and s2. If it finds the position, it returns the overlap length and the
 * overlap distance of the two strings. If not, it returns INT_MAX. If it finds two positions for both the
 * first and last window in s1, it chooses the one with the smaller distance.
 * Specifically, if the function finds the position, it returns three kinds of values.
 * 		If s2 is to the left of s1, the values(overlap length and distance) are negative integer;
 * 		If s2 is to the right of s1, the values(overlap length and distance) are positive integer;
 * 		If s1 and s2 has inclusion, the overlap distance is INT_MIN.
 *
 * @param s1 String the first string, s2 String the second string, d2Dis int d2 distance.
 * @return the first element is the overlap length, the second is the overlap distance.
 * If s2 is to the left of s1, the length are negative, the distance is zero or negative.
 * If s2 is to the right of s1, the length are positive, the distance is zero or positive.
 * If no overlap is found, the distance is INT_MAX.
 * If s2 is included in s1, the distance is INT_MIN.
 * If s1 is included in s2, the distance is INT_MIN.
 */
vector<int> OvlDistance::getOVLDistance(const string& tS1, const string& tS2, bool useHeuristic) {
	int s1Len = tS1.size();
	int s2Len = tS2.size();
	//We do not compute distances between sequences in non-adjacent groups
	int s1Flag = s1Len<=SHORT_LEN? 0 : (s1Len<=MEDIUM_LEN? 1 : 2);
	int s2Flag = s2Len<=SHORT_LEN? 0 : (s2Len<=MEDIUM_LEN? 1 : 2);
	vector<int> returnValues(2);
	if (abs(s1Flag-s2Flag)>1) {
		returnValues[0] = 0;
		returnValues[1] = INT_MAX;
		return returnValues;
	}

	int shorterLen = s1Len>s2Len? s2Len : s1Len;
	if (shorterLen <= SHORT_LEN) { //use parameters for short
		InclusionThreshold = INCLUSION_THRESHOLD_S;
	} else if (shorterLen <= MEDIUM_LEN) { //use parameters for medium
		InclusionThreshold = INCLUSION_THRESHOLD_M;
	} else { //use parameters for long
		InclusionThreshold = INCLUSION_THRESHOLD_L;
	}

	string s1 = "";
	string s2 = "";
	int flag = 1; //1 - no switch for tS1 and tS2; -1 - switch.
	/*
	 * put the shorter string to s1 and the longer one to s2 in
	 * order to identify inclusion. Now we just need to identify the
	 * situation when s1 is included in s2.
	 */
	if (tS1.length() > tS2.length()) {
		s1 = tS2;
		s2 = tS1;
		flag = -1; //tS1 and tS2 are switched
	} else {
		s1 = tS1;
		s2 = tS2;
	}

	BestWindowMatches best = d2.matchEndWindows(s1, s2, true, useHeuristic);
	if ((best.numBestLeftWindows==0) && (best.numBestRightWindows==0)) { //if not overlap, return
		returnValues[0] = 0;
		returnValues[1] = INT_MAX;
		return returnValues;
	}
	vector<int> tLeftPos = best.bestLeftStart;
	vector<int> tRightPos = best.bestRightStart;
	vector<int> leftPos = reducePos(tLeftPos, best.numBestLeftWindows);
	vector<int> rightPos = reducePos(tRightPos, best.numBestRightWindows);
	int windowSize = best.d2WindowSize;

	int disLeft = INT_MAX;
	int disRight = INT_MAX;
	int ovlDis = INT_MAX;
	int lenOverlap = 0;
	int lLenOverlap = 0;
	int rLenOverlap = 0;

	for (int i = 0; i < leftPos.size(); i++) {
		int lPos = leftPos[i];
		int tLenOverlap = s2.length() - lPos;
		int tmpDis = INT_MAX;
		if (tLenOverlap > s1.length()) { //if s1 is included in s2
			tmpDis = alignment.getDistance(s1.substr(0, s1.length()),
					s2.substr(lPos, s1.length()));
			tLenOverlap = s1.length();
		} else {
			tmpDis = alignment.getDistance(s1.substr(0, tLenOverlap),
					s2.substr(lPos));
		}
		if (tmpDis < disLeft) { // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
			disLeft = tmpDis;
			lLenOverlap = tLenOverlap;
		}
	}

	for (int i = 0; i < rightPos.size(); i++) {
		int rPos = rightPos[i];
		int tLenOverlap = rPos + windowSize;
		int lenInS1 = s1.length() - tLenOverlap;

		int tmpDis = INT_MAX;
		if (lenInS1 < 0) { //if s1 is included in s2
			tmpDis = alignment.getDistance(s1.substr(0), s2.substr(
					tLenOverlap - s1.length(), s1.length()));
			tLenOverlap = s1.length();
		} else {
			tmpDis = alignment.getDistance(s1.substr(lenInS1), s2.substr(
					0, tLenOverlap));
		}
		if (tmpDis < disRight) {// && (tLenOverlap > rLenOverlap). do we need to use two conditions or just one?
			disRight = tmpDis;
			rLenOverlap = tLenOverlap;
		}
	}

	// compare disLeft and disRight, select the one with smaller value.
	if (disLeft < disRight) {
		ovlDis = -1 * disLeft * flag; //minus represents that s2 is to the left of s1
		lenOverlap = -1 * lLenOverlap * flag;
	} else {
		ovlDis = disRight * flag; //s2 is to the right of s1
		lenOverlap = rLenOverlap * flag;
	}

	/*if s1 is included in s2, we will ignore s2 by assigning INT_MAX to the overlap distance.
	 * We do not need to consider that s1 includes s2 because we have switched them at the beginning
	 * of this function if s1 is longer than s2.
	 *
	 * Note: if ovlWindowSize > length of s1, disRight and disLeft would not be equal to zero.
	 */
	if ((abs(lenOverlap) == s1.length()) && (s1.length() <= s2.length())
			&& ((disRight <= InclusionThreshold) || (abs(disLeft)
					<= InclusionThreshold))) {
		lenOverlap = min(s1.length(), s2.length());
		ovlDis = INT_MIN;
	}

	returnValues[0] = lenOverlap;
	returnValues[1] = ovlDis;
	return returnValues;
}

vector<int> OvlDistance::reducePos(const vector<int>& input, int len) {
	vector<int> ret;
	if (len == 1) {
		ret.push_back(input[0]);
	} else if (len > 1) {
		ret.push_back(input[0]);
		ret.push_back(input[len - 1]);
	}
	return ret;
}
/*
 * judge if s1 is included in s2
 * @return true or false
 */
bool OvlDistance::checkInclusion(const string& s1, const string& s2, bool createNewHash) {
	int s1Len = s1.size();
	int s2Len = s2.size();
	//We do not compute distances between sequences in non-adjacent groups
	int s1Flag = s1Len<=SHORT_LEN? 0 : (s1Len<=MEDIUM_LEN? 1 : 2);
	int s2Flag = s2Len<=SHORT_LEN? 0 : (s2Len<=MEDIUM_LEN? 1 : 2);
	if (abs(s1Flag-s2Flag)>1) {
		return false;
	}

	int shorterLen = s1Len>s2Len? s2Len : s1Len;
	if (shorterLen <= SHORT_LEN) { //use parameters for short
		InclusionThreshold = INCLUSION_THRESHOLD_S;
	} else if (shorterLen <= MEDIUM_LEN) { //use parameters for medium
		InclusionThreshold = INCLUSION_THRESHOLD_M;
	} else { //use parameters for long
		InclusionThreshold = INCLUSION_THRESHOLD_L;
	}

	BestWindowMatches best = d2.matchEndWindows(s1, s2, createNewHash);
	if ((best.numBestLeftWindows==0) && (best.numBestRightWindows==0)) { //if not overlap, return false
		return false;
	}

	vector<int> tLeftPos = best.bestLeftStart;
	vector<int> tRightPos = best.bestRightStart;
	vector<int> leftPos = reducePos(tLeftPos, best.numBestLeftWindows);
	vector<int> rightPos = reducePos(tRightPos, best.numBestRightWindows);
	int windowSize = best.d2WindowSize;
	int disLeft = INT_MAX;
	int disRight = INT_MAX;
	int lenOverlap = 0;
	int lLenOverlap = 0;
	int rLenOverlap = 0;

	// if all leftPos[i] are -1, disLeft will be kept to be INT_MAX.
	for (int i = 0; i < leftPos.size(); i++) {
		int lPos = leftPos[i];
		int tLenOverlap = s2.length() - lPos;
		int tmpDis = INT_MAX;
		if (tLenOverlap > s1.length()) { //if s1 is included in s2
			tmpDis = alignment.getDistance(s1.substr(0, s1.length()),
					s2.substr(lPos, s1.length()));
			tLenOverlap = s1.length();
		} else {
			tmpDis = alignment.getDistance(s1.substr(0, tLenOverlap),
					s2.substr(lPos));
		}
		if (tmpDis < disLeft) { // && (tLenOverlap > lLenOverlap), do we need to use two conditions or just one?
			disLeft = tmpDis;
			lLenOverlap = tLenOverlap;
		}
	}

	// if all rightPos[i] are -1, disRight will be kept to be INT_MAX.
	for (int i = 0; i < rightPos.size(); i++) {
		int rPos = rightPos[i];
		int tLenOverlap = rPos + windowSize;
		int lenInS1 = s1.length() - tLenOverlap;

		int tmpDis = INT_MAX;
		if (lenInS1 < 0) { //if s1 is included in s2
			tmpDis = alignment.getDistance(s1.substr(0), s2.substr(tLenOverlap - s1.length(), s1.length()));
			tLenOverlap = s1.length();
		} else {
			tmpDis = alignment.getDistance(s1.substr(lenInS1), s2.substr(0, tLenOverlap));
		}
		if (tmpDis < disRight) {// && (tLenOverlap > rLenOverlap). do we need to use two conditions or just one?
			disRight = tmpDis;
			rLenOverlap = tLenOverlap;
		}
	}

	// compare disLeft and disRight, select the one with smaller value.
	if (disLeft < disRight) {
		lenOverlap = -1 * lLenOverlap;
	} else {
		lenOverlap = rLenOverlap;
	}

	/*if s1 is included in s2, return true, else return false
	 */
	if ((abs(lenOverlap) == s1.length()) && (s1.length() <= s2.length())
			&& ((disRight <= InclusionThreshold) || (abs(disLeft)
					<= InclusionThreshold))) {
		return true;
	} else {
		return false;
	}
}

