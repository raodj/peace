#ifndef ZY_Graph_HH
#define ZY_Graph_HH

#include <vector>
#include <string>
#include <limits.h>
#include <stack>
#include <map>
#include "Node.h"
#include "OvlDistance.h"
#include "CalculatedOvlDistance.h"
#include "InclusionNodes.h"
#include "SixTuple.h"
#include "Prim.h"
#include "Param.h"

class Graph {
public:
	int numOfLevels;
	int treeThreshold;
	std::vector<Node> graphNodes;
	OvlDistance ovl;
	CalculatedOvlDistance calDist;
	DefGraph mst; //minimum spanning tree generated from peace
	InclusionNodes* inc;

	Graph() {numOfLevels=0; inc=NULL;};
	Graph(InclusionNodes* in);
	std::vector<SixTuple*> get2CloseNodesFromMST();
	SixTuple get2CloseNodesFromGrand(int index, const SixTuple sixTuple);
	SixTuple checkLeftEndFromMST(int index, const SixTuple sixTuple);
	SixTuple checkRightEndFromMST(int index, const SixTuple sixTuple);
	std::vector<std::stack<int> > getNodesFromMST(std::vector<std::stack<int> > nodes); //also called by method in reconstruction.cpp
	inline void setMst(DefGraph& m) {
		mst = m;
	}
private:
	std::vector<SixTuple*> handleInclusion();
	SixTuple findAdjacentNode(std::stack<int> nodes, int index, const SixTuple sixTuple);

public:
	inline void addNode(Node s) {
		graphNodes.push_back(s);
	}

	inline void removeNode(int index) {
		graphNodes.erase(graphNodes.begin() + index);
	}

	inline int getSizeofGraph() {
		return graphNodes.size();
	}

	/*
	 * get length of the node with index i
	 */
	inline int getLenOfNode(int i) {
		return graphNodes.at(i).getLen();
	}

	/*
	 * get ID of the node with index i
	 */
	inline std::string getNameOfNode(int i) {
		return graphNodes.at(i).getName();
	}

	/*
	 * get sequence of the node with index i
	 */
	inline std::string getSeqOfNode(int i) {
		return graphNodes.at(i).getSeq();
	}

	/*
	 * get comment of the node with index i
	 */
	inline std::string getCommentOfNode(int i) {
		return graphNodes.at(i).getComment();
	}

	/*
	 * get quality scores of the node with index i
	 */
	inline std::vector<int> getQualScoresOfNode(int i) {
		return graphNodes.at(i).getQualScores();
	}

	/*
	 * set direction of the node.
	 * -1-complemented, 1-uncomplemented.
	 * only used for ACE format.
	 */
	inline void setDirectionOfNode(int nodeIdx, int direction) {
		graphNodes.at(nodeIdx).setDirection(direction);
	}

	inline std::string getDirectionOfNode(int nodeIdx) {
		return graphNodes.at(nodeIdx).direction;
	}
};

#endif
