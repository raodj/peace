#ifndef ZY_ESTAssembly_HH
#define ZY_ESTAssembly_HH

#include <ctype.h>
#include "Reconstruction.h"
#include "Graph.h"
#include "SixTuplesGeneration.h"
#include "InclusionNodes.h"
#include "Prim.h"

class ESTAssembly {
public:
	Graph* g;	//graph to store all the ests. It is generated in 'readEstFile' function.
	SixTuplesGeneration* gen;
	Reconstruction* rec;
	InclusionNodes* incNodes;

	ESTAssembly() {};
	ESTAssembly(std::string& estF, std::string& mstF);
	ESTAssembly(std::string& estF, std::string& mstF, std::string& qualF);
	void assemble(std::string& con, std::string& sing, std::string& numF);
	~ESTAssembly();

private:
	std::string estFileName;
	std::string mstFileName;
	std::string qualFileName;
	int readEstFile();
	int readEstQualFile();
	void readMST();
	std::vector<std::string> split(const std::string& str, char delimit);
	std::string toUpperCase(const std::string& str);
	std::vector<int> getIntScores(const std::string& str);
};

#endif
