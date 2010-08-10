#ifndef ZY_PARAM_HH
#define ZY_PARAM_HH

extern int SHORT_LEN;
extern int MEDIUM_LEN;
extern int WINDONW_SIZE_S;
extern int WINDONW_SIZE_M;
extern int WINDONW_SIZE_L;
extern int D2_THRESHOLD_S;
extern int D2_THRESHOLD_M;
extern int D2_THRESHOLD_L;
extern int D2_WORD_SIZE;
extern int U_S;
extern int U_M;
extern int U_L;
extern int T_S;
extern int T_M;
extern int T_L;
extern int UV_SKIP_S;
extern int UV_SKIP_M;
extern int UV_SKIP_L;
extern int TV_MAX_S;
extern int TV_MAX_M;
extern int TV_MAX_L;
extern int HEURISTIC_WORD_SIZE;
extern int INCLUSION_THRESHOLD_S;
extern int INCLUSION_THRESHOLD_M;
extern int INCLUSION_THRESHOLD_L;
extern int ALIGNMENT_THRESHOLD;
extern int NUMOFLEVELS;
extern int BAND_WIDTH_NW;
extern int USE_BOUNDED_NW;
extern int OUTPUT_SNP;
extern int USE_QUALITY_FILE;
extern int OUTPUT_ACE;
extern int LEN_DIFFERENCE;
extern int SHORTER_EST_LEN;

//the width of band which is used to speed up method "getSWAlignment".
//It should be <= COMPARISON_LENGTH. It is related to coverage depth, the bigger the coverage, the smaller this value.
const int BAND_WIDTH_SW=100;

//If use bounded version of Smith-Waterman algorithm to get SW Alignment when doing reconstruction.
//1-use, 0-not use. Use bounded version may reduce the precision.
//Bounded SW is not supported in current version.
const int USE_BOUNDED_SW=0;

#endif
