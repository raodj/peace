#ifndef ZY_InclusionNodes_HH
#define ZY_InclusionNodes_HH

#include <limits.h>
#include <string>
#include <vector>
#include <map>
#include "Prim.h"

class InclusionNodes {
public:
	std::map<int, int> nodes;
	std::multimap<int, int> nodes2;

	void addNode1(int chd, int parent);
	void addNode2(int chd, int parent);
	int getSize();
	bool containInclusionNode(int idx);
	std::vector<int> containPNode(int pIdx, int size);
	std::vector<int> getAllChdNodes();
	void printAllNodes();
private:
	DefGraph nodesGraph;
	void makeGraphFromAllNodes(int size);
};

#endif
