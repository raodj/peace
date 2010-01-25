#include <iostream>
#include <string>
#include "ESTAssembly.h"
#include "Param.h"
using namespace std;

int main(int argc, char* argv[]) {
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

	if (argc > 1) {
		ESTFILE = argv[1];
		MSTFILE = argv[2];
		CONSENSUSFILE = argv[3];
		SINGLETONFILE = argv[4];
		NUMOFUSEDEST = argv[5];
	}

	ESTAssembly assemble(ESTFILE, MSTFILE);
	assemble.assemble(CONSENSUSFILE, SINGLETONFILE, NUMOFUSEDEST);
	int nCall = assemble.rec->alignment.numCall;
	int uTime = assemble.rec->alignment.usedTime;
	int nCall2 = assemble.g->ovl.alignment.numCall2;
	int uTime2 = assemble.g->ovl.alignment.usedTime2;
	cout << "number of calls for SW alignment: " << nCall << endl;
	cout << "used time for SW alignment: " << uTime << endl;
	cout << "number of calls for NW score: " << nCall2 << endl;
	cout << "used time for NW score: " << uTime2 << endl;
}



