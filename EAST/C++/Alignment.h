#ifndef ZY_Alignment_HH
#define ZY_Alignment_HH

#include <limits.h>
#include <string>
#include <vector>
#include "AlignmentAlgorithm.h"
#include "Param.h"

class Alignment {
public:
	Alignment();
	~Alignment();

	int alignmentThreshold; //It's the threshold for alignment. That is, all the alignment with
							// the distance which is bigger than the value will be seen as infinity.
	int numCall;
	time_t usedTime;
	int numCall2;
	time_t usedTime2;
	int getDistance(const std::string& s1, const std::string& s2, std::vector<int> qualScore1=std::vector<int>(), std::vector<int> qualScore2=std::vector<int>());
	AlignResult getLocalAlignment(const std::string& s1, const std::string& s2, std::vector<int> qualScore1=std::vector<int>(), std::vector<int> qualScore2=std::vector<int>());
	AlignResult getBoundedLocalAlignment(const std::string& s1, const std::string& s2);
	void setScoringSystem(int match, int mismatch, int gap);

private:
	int getSimlarityScore(const std::string& s1, const std::string& s2);
	int getSimlarityScoreWithQual(const std::string& s1, const std::string& s2, const std::vector<int>& qualScore1, const std::vector<int>& qualScore2);
	AlignmentAlgorithm* alignAlgo;

};


#endif
