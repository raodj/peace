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
	intMatrix scoreMatrix; //A, C, G, T, N
	int gapPenalty;
	inline int encodeBase(char c) {
	    switch (c) {
	       case 'A' : return 0;
	       case 'C' : return 1;
	       case 'G' : return 2;
	       case 'T' : return 3;
	       case 'N' : return 4;
	       case 'P' : return 5;
	   }
	   return -1;
	}


public:
	AlignmentAlgorithm() {scoreMatrix = NULL;};
	AlignmentAlgorithm(int match, int mismatch, int gap);
	~AlignmentAlgorithm();
	inline void setScoringMatrix(const intMatrix& m) {scoreMatrix = m;}
	inline void setGapPenalty(int gap) {gapPenalty = gap;}
	int getNWScore(const std::string& s1, const std::string& s2);
	AlignResult getNWAlignment(const std::string& s1, const std::string& s2);
	int getSWScore(const std::string& s1, const std::string& s2);
	AlignResult getSWAlignment(const std::string& s1, const std::string& s2);
	AlignResult getBoundedSWAlignment(const std::string& horizontal, const std::string& vertical);
	std::vector<int> encodeString(const std::string& s1);
};



#endif
