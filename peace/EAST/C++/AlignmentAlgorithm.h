#ifndef ZY_AlignmentAlgorithm_HH
#define ZY_AlignmentAlgorithm_HH

#include <string>
#include <vector>
#include "Param.h"

typedef int** intMatrix;
intMatrix createIntMatrix(int n, int m);
void deleteIntMatrix(int** p, int n);

struct AlignResult {
	int score;
	std::string str1;
	std::string str2;
};

/*
 * This class implements Global and local alignment algorithm
 */
class AlignmentAlgorithm {
private:
	intMatrix scoreMatrix; //A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X
	int gapPenalty;
	inline int encodeBase(char c) {
	    switch (c) {
	       case 'A' : return 0;
	       case 'C' : return 1;
	       case 'G' : return 2;
	       case 'T' : return 3;
	       case 'N' : return 4;
	       case 'P' : return 5;
	       case 'Y' : return 6;
	       case 'R' : return 7;
	       case 'W' : return 8;
	       case 'S' : return 9;
	       case 'K' : return 10;
	       case 'M' : return 11;
	       case 'D' : return 12;
	       case 'V' : return 13;
	       case 'H' : return 14;
	       case 'B' : return 15;
	       case 'X' : return 16;
	   }
	   return -1;
	}


public:
	AlignmentAlgorithm() {scoreMatrix = NULL;};
	AlignmentAlgorithm(int match, int mismatch, int gap);
	~AlignmentAlgorithm();
	void setScoringMatrix(int match, int mismatch);
	inline void setGapPenalty(int gap) {gapPenalty = gap;}
	int getNWScore(const std::string& s1, const std::string& s2);
	int getBoundedNWScore(const std::string& s1, const std::string& s2);
	AlignResult getNWAlignment(const std::string& s1, const std::string& s2);
	int getNWScoreWithQual(const std::string& s1, const std::string& s2, const std::vector<int>& qualScore1, const std::vector<int>& qualScore2);
	int getBoundedNWScoreWithQual(const std::string& s1, const std::string& s2, const std::vector<int>& qual1, const std::vector<int>& qual2);

	int getSWScore(const std::string& s1, const std::string& s2);
	AlignResult getSWAlignment(const std::string& s1, const std::string& s2);
	AlignResult getSWAlignmentWithQual(const std::string& s1, const std::string& s2, const std::vector<int>& qualScore1, const std::vector<int>& qualScore2);
	AlignResult getBoundedSWAlignment(const std::string& s1, const std::string& s2);
	std::vector<int> encodeString(const std::string& s1);
};



#endif
