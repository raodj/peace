#include "SingleBase.h"
#include "Param.h"
#include <iostream>
using namespace std;

SingleBase::SingleBase(char c, int qualVal) {
	curBase = '\0';
	curBaseNum = 0;
	curqualVal = 0;
	numA = 0;
	numC = 0;
	numG = 0;
	numT = 0;
	numN = 0;
	numDash = 0;
	qualValA = 0;
	qualValC = 0;
	qualValG = 0;
	qualValT = 0;
	qualValN = 0;
	qualValDash = 0;
	changeBaseNum(c, qualVal);
	changeCurBase();
}

SingleBase::SingleBase(char c1, char c2, int qualVal1, int qualVal2) {
	curBase = '\0';
	curBaseNum = 0;
	curqualVal = 0;
	numA = 0;
	numC = 0;
	numG = 0;
	numT = 0;
	numN = 0;
	numDash = 0;
	qualValA = 0;
	qualValC = 0;
	qualValG = 0;
	qualValT = 0;
	qualValN = 0;
	qualValDash = 0;
	changeBaseNum(c1, qualVal1);
	changeBaseNum(c2, qualVal2);
	changeCurBase();
}

void SingleBase::addOneBase(char base, int qualVal) {
	changeBaseNum(base, qualVal);
	changeCurBase();
}

int SingleBase::getTotalNumOfBase() {
	return numA+numC+numG+numT+numN;
}

char SingleBase::getCurBase() {
	if (curBase == '-') {
		return 'P';
	} else {
		return curBase;
	}
}

void SingleBase::changeBaseNum(char base, int qualVal) {
	if (USE_QUALITY_FILE == 1) { //use quality file
		switch (base) {
		case 'A':
			numA++;
			qualValA += qualVal;
			break;
		case 'T':
			numT++;
			qualValT += qualVal;
			break;
		case 'C':
			numC++;
			qualValC += qualVal;
			break;
		case 'G':
			numG++;
			qualValG += qualVal;
			break;
		case 'N':
			numN++;
			qualValN += qualVal;
			break;
		case '-':
			numDash++;
			qualValDash += qualVal;
			break;
		}
	} else {
		switch (base) {
		case 'A':
			numA++;
			break;
		case 'T':
			numT++;
			break;
		case 'C':
			numC++;
			break;
		case 'G':
			numG++;
			break;
		case 'N':
			numN++;
			break;
		case '-':
			numDash++;
			break;
		}

	}
}

void SingleBase::changeCurBase() {
	if (USE_QUALITY_FILE == 1) { //use quality file
		int maxVal = qualValDash;
		int baseNum = numDash;
		char base = '-';
		if (maxVal < qualValC) {
			base = 'C';
			baseNum = numC;
			maxVal = qualValC;
		}
		if (maxVal < qualValG) {
			base = 'G';
			baseNum = numG;
			maxVal = qualValG;
		}
		if (maxVal < qualValA) {
			base = 'A';
			baseNum = numA;
			maxVal = qualValA;
		}
		if (maxVal < qualValT) {
			base = 'T';
			baseNum = numT;
			maxVal = qualValT;
		}
		if (maxVal < qualValN) {
			base = 'N';
			baseNum = numN;
			maxVal = qualValN;
		}
		curBase = base;
		curBaseNum = baseNum;
		curqualVal = maxVal;
	} else {
		int maxNum = numDash;
		char base = '-';

		//if the number of C/G is equal to A/T, we will take C/G.
		if (maxNum < numC) {
			base = 'C';
			maxNum = numC;
		}
		if (maxNum < numG) {
			base = 'G';
			maxNum = numG;
		}
		if (maxNum < numA) {
			base = 'A';
			maxNum = numA;
		}
		if (maxNum < numT) {
			base = 'T';
			maxNum = numT;
		}
		if (maxNum < numN) {
			base = 'N';
			maxNum = numN;
		}

		curBase = base;
		curBaseNum = maxNum;
	}
}


