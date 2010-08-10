#include <string>
#include "SepReverseCluster.h"
using namespace std;

int main(int argc, char* argv[]) {
	//The input ests are put into this file
	string ESTFILE = "estFile.fa";
	//The MST is in the file
	string MSTFILE = "mstFile.mst";

	if (argc > 1) {
		ESTFILE = argv[1];
		MSTFILE = argv[2];
	}


	SepReverseCluster sc(ESTFILE, MSTFILE);
	sc.operateC();
}



