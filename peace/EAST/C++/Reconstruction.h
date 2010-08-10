#ifndef ZY_Reconstruction_HH
#define ZY_Reconstruction_HH
// reconstruction from all the ESTs

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include <map>
#include <algorithm>
#include <stdlib.h>
#include "SixTuple.h"
#include "Prim.h"
#include "Alignment.h"
#include "InclusionNodes.h"
#include "Graph.h"
#include "SingleBase.h"
#include "Param.h"

class UsedNode {
public:
	int index; //index of the node
	int pos;
	UsedNode(int idx, int p);
	friend bool operator< (const UsedNode& k1, const UsedNode& k2);
};

class StartPos {
public:
	int pos;
	int index; //index of the node
	StartPos() {};
	StartPos(int p, int idx);
	friend bool operator< (const StartPos& k1, const StartPos& k2);
};

class LeftEnd {
public:
	int index;
	int lenOfSeq;
	std::string seq;
	LeftEnd() {};
	LeftEnd(int idx, std::string s);
	friend bool operator< (const LeftEnd& k1, const LeftEnd& k2);
};

class Reconstruction {
public:
	Graph* g;
	std::vector<SixTuple*> alignArray;
	std::vector<int> sPos;	//starting positions of all the nodes
				//the index in the array is the index of the node, the value is its starting position.
				//initialized in reconstructFromEnds method. That is, it is reset each time when the method is invoked.
	std::vector<int> sPosDebug;	//starting positions of all the nodes
						//it's used for debugging. All the left ends will be assigned to their actual value in order to calculate inversions later.
	Alignment alignment;
	InclusionNodes* incNodes;
	std::vector<SixTuple*> leftMostNodes;
	std::string consensusFileName;
	std::string singletonFileName;
	std::string numOfUsedESTsFileName;
	std::string printStr;
	std::vector<std::string> firstEsts;	//includes all the left end sequence corresponding to each element in allConsensus.
	std::vector<int> usedNodes;	//index of the array is the index of the node, 1-used, 0-not used. It is used to identify singletons.

	Reconstruction();
	Reconstruction(Graph* graph, std::vector<SixTuple*> align, std::vector<SixTuple*> leftEnds, InclusionNodes* inc, const std::string& con, const std::string& sing, const std::string& numF, int longestEstLen);
	void getConsensus();

private:
	int COMPARISON_LENGTH;
	int idxOfSNPFile; //used to record the index of the current output file for SNP analysis
	int idxOfContig; //used to record the index of the contig for ACE output
	int totalNumOfRD; //used to record the total number of reads in the ACE file
	void printConsensus();
	std::vector<std::string> reconstruct();
	std::vector<std::string> reconstuctForEachCluster();
	std::vector<std::vector<int> > genDGraph();
	std::vector<std::string> processLeftEnds(std::vector<int>& curLeftEnds, std::vector<std::vector<int> >& dGraph);
	std::string reconstructFromEnds(std::vector<std::vector<int> > leftEnds, std::vector<std::vector<int> >& dGraph);
	std::string processLeftEndsWithInclusion(std::vector<LeftEnd>& includeStrs, std::vector<std::vector<int> >& dGraph);
	int getNumUsedNodes(std::vector<StartPos>& a);
	std::string getCurConsensus(std::vector<SingleBase*> bases);
	std::vector<UsedNode> addInclusionNodes(std::vector<StartPos>& input);
	DefGraph constructMinTree(int source); //construct MST from the directed graph graphForPrim
	void getStartPos(int parentNode, DefGraph& tree, std::vector<std::vector<int> >& d);
	std::string replace(std::string str, const std::string& old, const std::string& newstr);
	std::string printLeftEndInfo(int leftEnd);
	std::vector<std::string> reconstructSeq(std::vector<StartPos>& a);
	std::vector<std::string> reconstructSeqWithQual(std::vector<StartPos>& a);
	void writeToSNPFile(std::vector<SingleBase*> bases, int totalLen);

	DefGraph graphForPrim; //this is a graph which is used for generating MST. It is generated once in genDGraph(), then used multiple times.
};



#endif
