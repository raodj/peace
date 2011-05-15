#ifndef ZY_PARAMDECLARE_HH
#define ZY_PARAMDECLARE_HH

#include <string>
/* Three groups of ESTs: short, medium, long
	short: length up to 150;
	medium: length 150-400;
	long: length 400 and up.
*/
int SHORT_LEN;
int MEDIUM_LEN;

//Parameters for D2

// The window size for d2
int WINDONW_SIZE_S;
int WINDONW_SIZE_M;
int WINDONW_SIZE_L;
//Threshold for d2. if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
//	THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
int D2_THRESHOLD_S;
int D2_THRESHOLD_M;
int D2_THRESHOLD_L;
//Use this value to define overlap distance of two inclusion subsequence. When there is no error in est, we can set it to be zero;
//  When error occurs, if the average overlap length is len, we can set it to be (1-(len-8)/len)*100; here 8 means  allowing 4 different
//   bases in two ests. That is, if the different bases <= 4, we assume them to be inclusion.
int INCLUSION_THRESHOLD_S;
int INCLUSION_THRESHOLD_M;
int INCLUSION_THRESHOLD_L;
//Bound of word for d2 distance. This is the word length we use in d2.
int D2_WORD_SIZE;
//The threshold for alignment. That is, all the alignment with the distance which is bigger than the value will be seen as infinity.
int ALIGNMENT_THRESHOLD;


//Parameters for D2 Heuristics
int HEURISTIC_WORD_SIZE;
int U_S;
int U_M;
int U_L;
int T_S;
int T_M;
int T_L;
int UV_SKIP_S;
int UV_SKIP_M;
int UV_SKIP_L;
int TV_MAX_S;
int TV_MAX_M;
int TV_MAX_L;

//The number of levels we will do in order to verify a left end. "0" means doing until leaves.
int NUMOFLEVELS;

//the width of band which is used to speed up method "getNWScore".
int BAND_WIDTH_NW;

//If use bounded version of Needleman-Wunsch algorithm to get NW score when getting overlap distance.
//1-use, 0-not use. Use bounded version may reduce the precision.
int USE_BOUNDED_NW;

//If generate output files for SNP analysis. 0: not, 1: generate.
int OUTPUT_SNP;
//The prefix of the output files for SNP analysis. The number of SNP output files is equal to the number of consensus sequences.
//They are named as snaAnalysis.0, snaAnalys.1, ...
//std::string SNP_FILE="snpAnalysis";

//If use quality files.
//1-use, 0-not use.
int USE_QUALITY_FILE;

//If output ACE file.
//1-yes, 0-no. If it's set 1, EAST will look for a file named as "$estFileName.aceinfo" and read into it.
int OUTPUT_ACE;

//This group of parameters are used for hybrid data.
//we will see the overlap distance of two ESTs to be infinite if
//	1)length difference > LEN_DIFFERENCE
//	2) and the shorter length < SHORTER_EST_LEN
int LEN_DIFFERENCE;
int SHORTER_EST_LEN;

// The progress file to which progress information is to be written
std::string PROGRESS_FILE_NAME;

#endif
