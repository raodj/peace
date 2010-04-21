#ifndef ZY_PARAM_HH
#define ZY_PARAM_HH

#include <string>
/* Three groups of ESTs: short, medium, long
	short: length up to 150;
	medium: length 150-400;
	long: length 400 and up.
*/
const int SHORT_LEN=150;
const int MEDIUM_LEN=400;

//Parameters for D2

// The window size for d2
const int WINDONW_SIZE_S=50;
const int WINDONW_SIZE_M=75;
const int WINDONW_SIZE_L=100;
//Threshold for d2. if the d2 distance is bigger than the threshold, we consider it to be infinite(=INT_MAX).
//	THRESHOLD = [(windowSize)-(boundOfWord)+1]^2
const int D2_THRESHOLD_S=2025;
const int D2_THRESHOLD_M=4900;
const int D2_THRESHOLD_L=9025;
//Use this value to define overlap distance of two inclusion subsequence. When there is no error in est, we can set it to be zero;
//  When error occurs, if the average overlap length is len, we can set it to be (1-(len-8)/len)*100; here 8 means  allowing 4 different
//   bases in two ests. That is, if the different bases <= 4, we assume them to be inclusion.
const int INCLUSION_THRESHOLD_S=27;
const int INCLUSION_THRESHOLD_M=27;
const int INCLUSION_THRESHOLD_L=22;
//Bound of word for d2 distance. This is the word length we use in d2.
const int D2_WORD_SIZE=6;
//The threshold for alignment. That is, all the alignment with the distance which is bigger than the value will be seen as infinity.
const int ALIGNMENT_THRESHOLD=40;


//Parameters for D2 Heuristics
const int HEURISTIC_WORD_SIZE=8;
const int U_S=4;
const int U_M=6;
const int U_L=8;
const int T_S=20;
const int T_M=30;
const int T_L=40;
const int UV_SKIP_S=4;
const int UV_SKIP_M=8;
const int UV_SKIP_L=8;
const int TV_MAX_S=50;
const int TV_MAX_M=75;
const int TV_MAX_L=100;

//The number of levels we will do in order to verify a left end. "0" means doing until leaves.
const int NUMOFLEVELS=0;

//The possible longest length of est. It's used in Reconstruction to reduce the length of strings doing local alignment
//const int COMPARISON_LENGTH=150;

//the width of band which is used to speed up method "getNWScore".
const int BAND_WIDTH_NW=20;

//If use bounded version of Needleman-Wunsch algorithm to get NW score when getting overlap distance.
//1-use, 0-not use. Use bounded version may reduce the precision.
const int USE_BOUNDED_NW=1;

//the width of band which is used to speed up method "getSWAlignment".
//It should be <= COMPARISON_LENGTH. It is related to coverage depth, the bigger the coverage, the smaller this value.
const int BAND_WIDTH_SW=100;

//If use bounded version of Smith-Waterman algorithm to get SW Alignment when doing reconstruction.
//1-use, 0-not use. Use bounded version may reduce the precision.
const int USE_BOUNDED_SW=0;

//If generate output files for SNP analysis. 0: not, 1: generate.
const int OUTPUT_SNP=0;
//The prefix of the output files for SNP analysis. The number of SNP output files is equal to the number of consensus sequences.
//They are named as snaAnalysis.0, snaAnalys.1, ...
//std::string SNP_FILE="snpAnalysis";
#endif
