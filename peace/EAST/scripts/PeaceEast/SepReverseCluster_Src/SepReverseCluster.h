#ifndef ZY_SepReverseCluster_HH
#define ZY_SepReverseCluster_HH

#include <ctype.h>
#include <vector>
#include <string>
#include <stack>
#include <iostream>
#include <fstream>
#include "Node.h"
#include "Prim.h"

class SepReverseCluster {
public:
	std::vector<Node> graphNodes;
	DefGraph mst; //minimum spanning tree generated from peace

	SepReverseCluster(const std::string& estF, const std::string& mstF);
	void operateC();
	void writeACluster(std::vector<int>& nodesInCluster, const std::string& newEstFile); //extract the EST file for the cluster that includes the node idx, then we can use peace to generate a new MST file.

private:
	std::string estFileName;
	std::string mstFileName;
	std::string newEstFileName;
	float treeThreshold;
	std::ofstream aceOutFile;
	void readEstFile(const std::string& inFileName);
	void readMST(const std::string& inFileName);
	void ProcessEachCluster();
	std::vector<std::stack<int> > getNodesFromMST(std::vector<std::stack<int> > nodes);
	std::vector<std::string> split(const std::string& str, char delimit);
	std::string toUpperCase(const std::string& str);
};

#endif
