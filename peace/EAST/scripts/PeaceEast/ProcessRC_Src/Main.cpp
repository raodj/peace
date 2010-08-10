#include <string>
#include "ProcessRC.h"
using namespace std;

int main(int argc, char* argv[]) {
	//The input ests are put into this file
	string ESTFILE = "estFile.fa";
	//The MST is in the file
	string MSTFILE = "mstFile.mst";
	//The new EST file
	string NEWMSTFILE = "NewEstFile.fa";

	if (argc > 1) {
		ESTFILE = argv[1];
		MSTFILE = argv[2];
		NEWMSTFILE = argv[3];
	}


	ProcessRC rc(ESTFILE, MSTFILE, NEWMSTFILE);
	rc.operateRC();
}



