#include <iostream>
#include <string>
#include "ESTAssembly.h"
#include "AlignmentAlgorithm.h"
#include "ParamDeclare.h"
#include "ProgressReporter.h"

using namespace std;

void setParamDefault() {
	SHORT_LEN=150;
	MEDIUM_LEN=400;
	WINDONW_SIZE_S=50;
	WINDONW_SIZE_M=75;
	WINDONW_SIZE_L=100;
	D2_THRESHOLD_S=2025;
	D2_THRESHOLD_M=4900;
	D2_THRESHOLD_L=9025;
	INCLUSION_THRESHOLD_S=27;
	INCLUSION_THRESHOLD_M=27;
	INCLUSION_THRESHOLD_L=22;
	D2_WORD_SIZE=6;
	ALIGNMENT_THRESHOLD=40;
	HEURISTIC_WORD_SIZE=8;
	U_S=4;
	U_M=6;
	U_L=8;
	T_S=20;
	T_M=30;
	T_L=40;
	UV_SKIP_S=4;
	UV_SKIP_M=8;
	UV_SKIP_L=8;
	TV_MAX_S=50;
	TV_MAX_M=75;
	TV_MAX_L=100;
	NUMOFLEVELS=0;
	BAND_WIDTH_NW=20;
	USE_BOUNDED_NW=1;
	OUTPUT_SNP=0;
	USE_QUALITY_FILE=0;
	OUTPUT_ACE=0;
	LEN_DIFFERENCE=300;
	SHORTER_EST_LEN=300;
    PROGRESS_FILE_NAME="";
}

int parse_args(int argc, char* argv[]) {
	int marker = 1;    // Remember that argv[0] is the executable name, so we don’t want to look there.

	while (marker < argc && argv[marker][0] == '-') {   // We assume an optional switch always starts with a “-’
		string switchStr = argv[marker];
		if (switchStr == "-SHORT_LEN") {
			SHORT_LEN = atoi(argv[++marker]);
		} else if (switchStr == "-MEDIUM_LEN") {
			MEDIUM_LEN = atoi(argv[++marker]);
		} else if (switchStr == "-WINDONW_SIZE_S") {
            WINDONW_SIZE_S = atoi(argv[++marker]);
        } else if (switchStr == "-WINDONW_SIZE_M") {
            WINDONW_SIZE_M = atoi(argv[++marker]);
        } else if (switchStr == "-WINDONW_SIZE_L") {
            WINDONW_SIZE_L = atoi(argv[++marker]);
        } else if (switchStr == "-D2_THRESHOLD_S") {
			D2_THRESHOLD_S = atoi(argv[++marker]);
		} else if (switchStr == "-D2_THRESHOLD_M") {
			D2_THRESHOLD_M = atoi(argv[++marker]);
		} else if (switchStr == "-D2_THRESHOLD_L") {
			D2_THRESHOLD_L = atoi(argv[++marker]);
		} else if (switchStr == "-INCLUSION_THRESHOLD_S") {
			INCLUSION_THRESHOLD_S = atoi(argv[++marker]);
		} else if (switchStr == "-INCLUSION_THRESHOLD_M") {
			INCLUSION_THRESHOLD_M = atoi(argv[++marker]);
		} else if (switchStr == "-INCLUSION_THRESHOLD_L") {
			INCLUSION_THRESHOLD_L = atoi(argv[++marker]);
		} else if (switchStr == "-D2_WORD_SIZE") {
			D2_WORD_SIZE = atoi(argv[++marker]);
		} else if (switchStr == "-ALIGNMENT_THRESHOLD") {
			ALIGNMENT_THRESHOLD = atoi(argv[++marker]);
		} else if (switchStr == "-HEURISTIC_WORD_SIZE") {
			HEURISTIC_WORD_SIZE = atoi(argv[++marker]);
		} else if (switchStr == "-U_S") {
			U_S = atoi(argv[++marker]);
		} else if (switchStr == "-U_M") {
			U_M = atoi(argv[++marker]);
		} else if (switchStr == "-U_L") {
			U_L = atoi(argv[++marker]);
		} else if (switchStr == "-T_S") {
			T_S = atoi(argv[++marker]);
		} else if (switchStr == "-T_M") {
			T_M = atoi(argv[++marker]);
		} else if (switchStr == "-T_L") {
			T_L = atoi(argv[++marker]);
		} else if (switchStr == "-UV_SKIP_S") {
			UV_SKIP_S = atoi(argv[++marker]);
		} else if (switchStr == "-UV_SKIP_M") {
			UV_SKIP_M = atoi(argv[++marker]);
		} else if (switchStr == "-UV_SKIP_L") {
			UV_SKIP_L = atoi(argv[++marker]);
		} else if (switchStr == "-TV_MAX_S") {
			TV_MAX_S = atoi(argv[++marker]);
		} else if (switchStr == "-TV_MAX_M") {
			TV_MAX_M = atoi(argv[++marker]);
		} else if (switchStr == "-TV_MAX_L") {
			TV_MAX_L = atoi(argv[++marker]);
		} else if (switchStr == "-LEN_DIFFERENCE") {
			LEN_DIFFERENCE = atoi(argv[++marker]);
		} else if (switchStr == "-SHORTER_EST_LEN") {
			SHORTER_EST_LEN = atoi(argv[++marker]);
		} else if (switchStr == "-NUMOFLEVELS") {
			NUMOFLEVELS = atoi(argv[++marker]);
		} else if (switchStr == "-USE_BOUNDED_NW") {
			USE_BOUNDED_NW = 1;
		} else if (switchStr == "-BAND_WIDTH_NW") {
			BAND_WIDTH_NW = atoi(argv[++marker]);
		} else if (switchStr == "-OUTPUT_SNP") {
			OUTPUT_SNP = 1;
		} else if (switchStr == "-q" || switchStr == "-USE_QUALITY_FILE") {
			USE_QUALITY_FILE = 1;
		} else if (switchStr == "-ace" || switchStr == "-OUTPUT_ACE") {
			OUTPUT_ACE = 1;
		} else if (switchStr == "-CONVERT_TO_SAM") {
            // Nothing to be done. We just consume this parameter
            // here. It is really used in the east.sh wrapper script
            // that is used by the GUI.
        } else if (switchStr == "--progress") {
            PROGRESS_FILE_NAME = argv[++marker];
        }
		marker++;
	}
	return marker;   // The first non-optional parameter will now be at argv[marker]
}


int main(int argc, char* argv[]) {
	setParamDefault(); //set default values for all parameters
	int marker = parse_args(argc, argv);

	//The input ests are put into this file
	string ESTFILE = "estFile.fa";
	//The MST is in the file
	string MSTFILE = "mstFile.mst";
	//Consensus sequencs are put into this file
	string CONSENSUSFILE = "consensus.out";
	//Singletons are put into this file
	string SINGLETONFILE = "singleton.out";
	//The number of ESTS that EAST used in the assembly is put into the file
	string NUMOFUSEDEST = "numOfUsedEsts.out";
	//Quality file
	string QUALITYFILE = "estFile.fa.qual";

	if (marker < argc) {
		ESTFILE = argv[marker];
		MSTFILE = argv[marker+1];
		CONSENSUSFILE = argv[marker+2];
		SINGLETONFILE = argv[marker+3];
		NUMOFUSEDEST = argv[marker+4];
		if (USE_QUALITY_FILE == 1) {
			QUALITYFILE = argv[marker+5];
		}
	}

    // Initialize progress reporting if requested by user
    ProgressReporter& pr = ProgressReporter::get();    
    if (!PROGRESS_FILE_NAME.empty()) {
        pr.initialize(PROGRESS_FILE_NAME, 4);
        pr.reportProgress(0);
    }
    
	ESTAssembly* assemble = NULL;

	if (USE_QUALITY_FILE == 1) { //use quality file
		assemble = new ESTAssembly(ESTFILE, MSTFILE, QUALITYFILE);
	} else {
		assemble = new ESTAssembly(ESTFILE, MSTFILE);
	}

	assemble->assemble(CONSENSUSFILE, SINGLETONFILE, NUMOFUSEDEST);
	int nCall = assemble->rec->alignment.numCall;
	int uTime = assemble->rec->alignment.usedTime;
	int nCall2 = assemble->g->ovl.alignment.numCall2;
	int uTime2 = assemble->g->ovl.alignment.usedTime2;
	cout << "number of calls for SW alignment: " << nCall << endl;
	cout << "used time for SW alignment: " << uTime << endl;
	cout << "number of calls for NW score: " << nCall2 << endl;
	cout << "used time for NW score: " << uTime2 << endl;

	delete assemble;
    // Report final step of progress
    pr.reportProgress(4);
    // All done successfully.
    return 0;
}
