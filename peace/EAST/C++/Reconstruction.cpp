#include "Reconstruction.h"
#include <sstream>
#include <cstdlib>

using namespace std;

UsedNode::UsedNode(int idx, int p) {
	index = idx;
	pos = p;
}

bool operator< (const UsedNode& k1, const UsedNode& k2) {
	if (k1.index < k2.index)
		return true;
	else if ((k1.index == k2.index) && (k1.pos < k2.pos))
		return true;
	else
		return false;
}

StartPos::StartPos(int p, int idx) {
	pos = p;
	index = idx;
}

bool operator< (const StartPos& k1, const StartPos& k2) {
	return (k1.pos < k2.pos);
}

LeftEnd::LeftEnd(int idx, string s) {
	index = idx;
	lenOfSeq = s.length();
	seq = s;
}

bool operator< (const LeftEnd& k1, const LeftEnd& k2) {
	return (k1.lenOfSeq < k2.lenOfSeq);
}

Reconstruction::Reconstruction() {
	g = NULL;
	incNodes = NULL;
}

Reconstruction::Reconstruction(Graph* graph, vector<SixTuple*> align, vector<SixTuple*> leftEnds, InclusionNodes* inc, const std::string& con, const std::string& sing, const std::string& numF) {
	g = graph;
	alignArray = align;
	incNodes = inc;
	leftMostNodes = leftEnds;
	consensusFileName = con;
	singletonFileName = sing;
	numOfUsedESTsFileName = numF;

	usedNodes = vector<int> (g->graphNodes.size(), 0);
}

void Reconstruction::getConsensus() {
	time_t start,end;
	cout << "Start to reconstruct." << endl;
	start = time(NULL);
	printConsensus();
	end = time(NULL);
	cout << "End to reconstruct." << endl;
	cout << "The time used to reconstruct is " << end-start << endl;
}

/*
 * Print the original sequence and multiple consensus into a file which is specified in the property file.
 */
void Reconstruction::printConsensus() {
	/*
	 * print consensus sequences
	 */
	ofstream outFile1;
	outFile1.open(consensusFileName.c_str(),ios::trunc);

	vector<string> consensus = reconstruct();
	int index = 1;
	for (int i=0; i<consensus.size(); i++) {
		string str = consensus[i];
		size_t found = str.find('\n');

		if (found != string::npos) { //there is "\n" in the sequence
			outFile1 << ">contig " << (index++) << endl;

			int start = 0;
			while (found != string::npos) {
				outFile1 << str.substr(start, (int)found-start) << endl;
				start = found + 1;
				found = str.find('\n', found+1);
			}
			outFile1 << str.substr(start) << endl;
		} else {
			outFile1 << ">contig " << (index++) << endl;
			outFile1 << consensus[i] << endl;
		}
	}
	outFile1.close();


	/*
	 * print singletons and numOfUsedESTs
	 */
	ofstream outFile2;
	outFile2.open(singletonFileName.c_str(),ios::trunc);

	int num = 0;
	for (int i=0; i<usedNodes.size(); i++) {
		if (usedNodes[i] == 0) {	//singleton
			outFile2 << ">" << g->getCommentOfNode(i) << "\n";
			outFile2 << g->getSeqOfNode(i) << "\n";
		} else {
			num++;
		}
	}
	outFile2.close();

	ofstream outFile3;
	outFile3.open(numOfUsedESTsFileName.c_str(),ios::trunc);
	outFile3 << ">number of used ESTs by EAST\n";
	outFile3 << num;
	outFile3.close();
}

/*
 * reconstruct the sequence
 *
 * @return an arraylist: The assembled sequences.
 */
std::vector<std::string> Reconstruction::reconstruct() {
	sPosDebug = vector<int> (g->graphNodes.size());

	vector<string> ret = processLeftEnds();

	//for debug
	std::stringstream out;
	out << printStr << ret.size() << " consensus from above " << leftMostNodes.size() << " left ends:\n";
	printStr = out.str();
	for (int p=0; p<ret.size(); p++) {
		printStr = printStr + ret[p] + "\n";
	}

	//print debug information about the generated consensus.
	cout << "*********************consensus:********************\n";
	cout << printStr << endl;

	return ret;
}

std::vector<std::vector<int> > Reconstruction::genDGraph() {
	/*
	 * Calculate the length of dGraph.
	 */
	int len = 0;
	for (int i=0; i<alignArray.size(); i++) {
		if (alignArray[i]->leftNode != -1) {
			len++;
		}
		if (alignArray[i]->rightNode != -1) {
			len++;
		}
	}

	/*
	 * generate dGraph.
	 */
	vector<vector<int> > dGraph = vector<vector<int> > (len, std::vector<int>(4));
	int indexOfDGraph = 0;
	for (int i=0; i<alignArray.size(); i++) {
		SixTuple* curTuple = alignArray[i];
		if (curTuple->leftNode != -1) {
			dGraph[indexOfDGraph][0] = curTuple->leftNode;
			dGraph[indexOfDGraph][1] = curTuple->curNode;
			dGraph[indexOfDGraph][2] = abs(curTuple->lDis);	//distance
			dGraph[indexOfDGraph][3] = abs(curTuple->lOvlLen);	//overlap length
			indexOfDGraph++;
		}

		if (curTuple->rightNode != -1) {
			dGraph[indexOfDGraph][0] = curTuple->curNode;
			dGraph[indexOfDGraph][1] = curTuple->rightNode;
			dGraph[indexOfDGraph][2] = curTuple->rDis;	//distance
			dGraph[indexOfDGraph][3] = curTuple->rOvlLen;	//overlap length
			indexOfDGraph++;
		}
	}

	return dGraph;
}



/*
 *  Return consensus from this set of left ends.
 *  input: leftEnds[][2] : [][0]-index, [][1]-position of the left end; size-the size of the array
 *
 * 1. Generate dGraph.
 * 2. For each left-end node, starting from it to calculate positions for each node.
 * In order to get starting positions, it constructs a MST. The weight of the Minimum Spanning tree is
 * the overlap distance instead of overlap length.
 * Then reconstruct the sequence from the set of ESTs.
 *
 */
string Reconstruction::reconstructFromEnds(vector<vector<int> > leftEnds, vector<vector<int> >& dGraph) {
	string ret = "";
	int size = leftEnds.size();
	map<int, int> ends;
	map<int, int>::iterator ite;
	sPos = vector<int> (g->graphNodes.size(), 0);	//store starting positions of all the nodes, reset it to the initial state

	for (int i=0; i<size; i++) {
		int leftEnd = leftEnds[i][0];
		sPos[leftEnd] = leftEnds[i][1];
		ends[leftEnd] = 0;
		// Calculate starting positions using minimum spanning tree starting from this left-end node.
		DefGraph primMST = constructMinTree(g->graphNodes.size(), dGraph, leftEnd); //the first param is the total number of ESTs.

		//get starting positions for the nodes in primMST
		getStartPos(leftEnd, primMST, dGraph);

		printStr = printStr + printLeftEndInfo(leftEnd); //get information of this left end which is used to start reconstruction.
	}
	stringstream tstr;
	tstr << leftEnds[0][0];
	printStr = printStr + "Start from node " + tstr.str() + " to do reconstruction" + "\n\n";


	vector<StartPos> tmpArray;
	for (int j=0; j<sPos.size(); j++) {
		ite = ends.find(j);
		if ((ite != ends.end()) || (sPos[j] != 0)) { // is left end or a node which has got its position
			tmpArray.push_back(StartPos(sPos[j], j));
		}
	}


	vector<string> tStr = reconstructSeq(tmpArray);
	if (tStr.size() != 0) {
		printStr = printStr + tStr[1] + " nodes are used to reconstruct the sequence.\n";
		printStr = printStr + tStr[0] + "\n\n";
		ret = tStr[0];
	}

	return ret;
}

/*
 * This method is designed to group left ends into independent sets. All the dependent left ends (include one by one) are put
 * into one group.
 * In each group, we find the left end with the longest length and the one which will use the most number of ESTs to reconstruct,
 * and combine them together to form a consensus by calling the function "processLeftEndsWithInclusion".

 * @return all the consensus sequences.
 */
vector<string>Reconstruction::processLeftEnds() {
	vector<string> allOutputContigs;	//store all the generated sequences

	vector<vector<int> > dGraph = genDGraph();
	int sizeOfs = leftMostNodes.size();
	if (sizeOfs == 0) {
		cout << "zero left ends!" << endl;
		exit(1);
	}
	vector<LeftEnd> resultArray(sizeOfs); //store the starting positions of ests
	for (int i=0; i<sizeOfs; i++) { //start for
		int idx = leftMostNodes[i]->curNode;
		resultArray[i] = LeftEnd(idx, g->getSeqOfNode(idx));
	} //end for
	sort(resultArray.begin(), resultArray.end());

	vector<LeftEnd> allLeftEnds;
	for (int i=sizeOfs-1; i>=0; i--) {
		allLeftEnds.push_back(resultArray[i]);
	}

	while (true) {
		string s1 = allLeftEnds[0].seq;
		int s1Idx = allLeftEnds[0].index;
		vector<LeftEnd> includedEnds;
		includedEnds.push_back(allLeftEnds[0]);

		vector<LeftEnd> excludedEnds;
		for (int i=1; i<allLeftEnds.size(); i++) {
			vector<int> ovlDis = g->calDist.searchDistance(allLeftEnds[i].index, s1Idx);

			if (ovlDis[1] == INT_MIN) { //if resultArray[i].firstEst is included in s1
				includedEnds.push_back(allLeftEnds[i]);
				incNodes->addNode2(allLeftEnds[i].index, s1Idx);
			} else if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) { //have not been calculated
				bool b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1);
				if (b) {
					includedEnds.push_back(allLeftEnds[i]);
					incNodes->addNode2(allLeftEnds[i].index, s1Idx);
				}
				else
					excludedEnds.push_back(allLeftEnds[i]);
			} else { //have been calculated, but no inclusion
				excludedEnds.push_back(allLeftEnds[i]);
			}
		}

		string contig = processLeftEndsWithInclusion(includedEnds, dGraph);
		if ((contig.size() != 0) && (contig.erase(contig.find_last_not_of(" \n\r\t")+1).size() != 0)) { //trim
			allOutputContigs.push_back(contig);
		}
		if (excludedEnds.size() == 0) {
			break;
		} else {
			allLeftEnds = excludedEnds;
		}
	}
	return allOutputContigs;
}

/*
 * This method is called by "processLeftEnds".
 * This method is used to process those left ends that include each other.
 *
 * We combine the one with the longest first EST and the one using the most number of ESTs together to
 * form a new one and return it.
 */
string Reconstruction::processLeftEndsWithInclusion(vector<LeftEnd>& includeStrs, vector<vector<int> >& dGraph) {
	if (includeStrs.size() == 0) {
		return "";
	} else if (includeStrs.size() == 1) {
		vector<vector<int> > leftEnds = vector<vector<int> > (1, std::vector<int>(2));
		leftEnds[0][0] = includeStrs[0].index; //index
		leftEnds[0][1] = 0; //position
		string s = reconstructFromEnds(leftEnds, dGraph);
		return s;
	} else {
		int maxLen = 0; //the maximal length of the first EST.
		int idxMaxLen = 0;
		int numOfEnds = includeStrs.size();
		vector<vector<int> > leftEnds = vector<vector<int> > (numOfEnds, std::vector<int>(2));

		for (int i=0; i<numOfEnds; i++) {
			int tLen = includeStrs[i].lenOfSeq;
			if (tLen > maxLen) {
				maxLen = tLen;
				idxMaxLen = i;
			}
		}

		int idx = 0;
		string s1 = includeStrs[idxMaxLen].seq;
		leftEnds[idx][0] = includeStrs[idxMaxLen].index; //index
		leftEnds[idx++][1] = 0; //position

		for (int i=0; i<numOfEnds; i++) {
			if (i == idxMaxLen) continue;
			string s2 = includeStrs[i].seq;
			AlignResult strs = alignment.getLocalAlignment(s1, s2);
			int offset1 = s1.find(replace(strs.str1, "-", ""));
			int offset2 = s2.find(replace(strs.str2, "-", ""));
			int endPos = (offset1-offset2) > 0 ? (offset1-offset2) : 0;
			leftEnds[idx][0] = includeStrs[i].index;
			leftEnds[idx++][1] = endPos;
		}

		return reconstructFromEnds(leftEnds, dGraph);
	}
}

/*
 * reconstruct a sequence which starts from a left end.
 * return: ret[0]-the consensus, ret[1]-the total number of nodes used for reconstruction.
 */
vector<string> Reconstruction::reconstructSeq(vector<StartPos>& a) {
	vector<string> ret(2);
	int sizeOfa = a.size();
	if ((sizeOfa == 0) || (sizeOfa == 1)){
		//size=1, this node should be a singleton if it is not used in other place. it shouldn't appear as a consensus sequence.
		//size=1, then usedNodes[this node] will not be set to 1 because addInclusionNodes won't be called.
		return ret;
	}

	sort(a.begin(), a.end());
	vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction and set usedNodes[related nodes]=1.
	std::stringstream out;
	out << addedNodes.size();
	ret[1] = out.str();

	vector<StartPos> resultArray(addedNodes.size());
	for (int i=0; i<addedNodes.size(); i++) {
		UsedNode tmpNode = addedNodes[i];
		resultArray[i] = StartPos(tmpNode.pos, tmpNode.index);
	}
	sort(resultArray.begin(), resultArray.end());

	int comparisonLen = COMPARISON_LENGTH;
	vector<SingleBase*> bases;
	string tConsensus = g->getSeqOfNode(resultArray[0].index);
	string curSeq = "";
	int len = resultArray.size() - 1;
	for (int i=1; i<=len; i++) {
		curSeq = g->getSeqOfNode(resultArray[i].index);
		string tmpConsensus = tConsensus;
		if (tmpConsensus.length() > comparisonLen) {
			tmpConsensus = tConsensus.substr(tConsensus.size()-comparisonLen+1);
		}

		AlignResult strs;
		if (USE_BOUNDED_SW == 0) { //use ordinary version
			strs = alignment.getLocalAlignment(tmpConsensus, curSeq);
		} else  { //use bounded version
			strs = alignment.getBoundedLocalAlignment(tmpConsensus, curSeq);
		}
		tConsensus = replace(tConsensus, replace(strs.str1, "-", ""), strs.str1);
		int offset = (int)tConsensus.find(strs.str1);

		string tSeq = replace(curSeq, replace(strs.str2, "-", ""), strs.str2);
		int tmpOff = (int)tSeq.find(strs.str2);
		string firstPartCur = tSeq.substr(0, tmpOff);
		curSeq = tSeq.substr(tmpOff);
		if (i == 1) {
			int len1 = tConsensus.size();
			int len2 = curSeq.size();
			int end = max(offset+len2, len1);
			for (int j=0; j<end; j++) {
				if ((j < len1) && (j-offset >= 0) && (j-offset < len2)) { //overlap part
					bases.push_back(new SingleBase(tConsensus[j], curSeq[j-offset]));
				} else if ((j-offset < 0) || (j-offset >= len2)) {
					bases.push_back(new SingleBase(tConsensus[j]));
				} else if (j >= len1) {
					bases.push_back(new SingleBase(curSeq[j-offset]));
				}
			}
		} else {
			int len1 = tConsensus.size();
			int len2 = curSeq.size();
			int end = max(offset+len2, len1);
			for (int j=offset; j<end; j++) {
				if ((j < len1) && (j-offset < len2)) { //overlap part
					char c1 = tConsensus[j];
					char c2 = curSeq[j-offset];

					if (c1 != '-') {
						bases[j]->addOneBase(c2);
					} else {
						bases.insert(bases.begin()+j, new SingleBase(c1, c2));
					}
				} else if (j >= len1) {
					bases.push_back(new SingleBase(curSeq[j-offset]));
				}
			}

			//put first part of curSeq into bases
			for (int k=0; k<tmpOff; k++) {
				if (offset-tmpOff+k >= 0)
					bases[offset-tmpOff+k]->addOneBase(firstPartCur[k]);
			}
		}

		tConsensus = getCurConsensus(bases);
	}

	ret[0] = replace(tConsensus, "P", "");

	//release bases
	for (int i=0; i<bases.size(); i++) {
		delete bases[i];
	}
	return ret;
}

int Reconstruction::getNumUsedNodes(vector<StartPos>& a) {
	int sizeOfa = a.size();
	if (sizeOfa == 0) {
		return 0;
	} else if (sizeOfa == 1) {
		return 1;
	}

	sort(a.begin(), a.end());

	vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction.
	return addedNodes.size();
}


string Reconstruction::getCurConsensus(vector<SingleBase*> bases) {
	int len = bases.size();
	string tStr;
	for (int i=0; i<len; i++) {
		tStr.push_back(bases[i]->getCurBase());
	}
	return tStr;
}

/*
 * in string 'str', replace 'old' with 'new'.
 */
string Reconstruction::replace(string str, const string& old, const string& newstr) {
	size_t found = str.find(old);
	while (found != string::npos) { //there is "old" in the sequence
		str.replace(found, old.size(), newstr);
		found = str.find(old, found+old.size()-1);
	}
	return str;
}


/*
 * Add all the inclusion nodes into the input arraylist.
 * For each element in the arraylist, put its corresponding node just after it.
 */
vector<UsedNode> Reconstruction::addInclusionNodes(vector<StartPos>& input) {
	map<int, int> tmpList;
	vector<UsedNode> retList;
	int size = input.size();

	for (int i=0; i<size; i++) {
		int curIdx = input[i].index;
		usedNodes[curIdx] = 1; //mark this node being used
		int pos = input[i].pos;
		tmpList.insert(pair<int, int> (curIdx, pos));

		vector<int> chdIdx = incNodes->containPNode(curIdx, g->graphNodes.size()); //inclusion children index of the curIdx if exist.
		for (int j=0; j<chdIdx.size(); j++) {
			tmpList.insert(pair<int, int> (chdIdx[j], pos+1));
			usedNodes[chdIdx[j]] = 1; //mark this node being used
		}
	}
	for (map<int, int>::iterator it = tmpList.begin(); it != tmpList.end(); it++) {
		retList.push_back(UsedNode(it->first, it->second));
	}
	sort(retList.begin(), retList.end());
	return retList;
}

/*
 * Construct a directed Miminum spanning tree.
 *
 *  @param nOfNodes number of nodes
 *  @param g a directed graph, the second dimension has three elements:
 *  	index of starting node, index of ending node, weight between them.
 */
DefGraph Reconstruction::constructMinTree(int nOfNodes, const vector<vector<int> >& input, int source) {
	// Make a directed graph.
	DefGraph dGraph(nOfNodes);
	for (int j=0; j<input.size(); j++) {
		if (input[j][3] != 0) {	//there is an edge between the nodes
			dGraph[input[j][0]].push_back(Edge(input[j][1], input[j][2]));
		}
	}
	DefGraph mst = Prim(dGraph, source, true);
	return mst;
}

/*
 * Calculate starting positions for each node.
 */
void Reconstruction::getStartPos(int parentNode, DefGraph& tree, vector<vector<int> >& d) {
	for (EdgeIterator it=tree[parentNode].begin(); it!=tree[parentNode].end(); it++) {
		int index = (*it).node;

		int overlapLen = 0;
		for (int i=0; i<d.size(); i++) {
			if ((d[i][0] == parentNode) && (d[i][1] == index)) {
				overlapLen = d[i][3];
				break;
			}
		}

		sPos[index] = sPos[parentNode] + g->getLenOfNode(parentNode) - overlapLen;
		getStartPos(index, tree, d);
	}
}


string Reconstruction::printLeftEndInfo(int leftEnd) {
	string ret = "";
	std::stringstream out;
	out << "Node " << leftEnd << " starts from " << g->getNameOfNode(leftEnd) << endl;
	ret = out.str();
	out.str("");
	out.clear();

	int sp = atoi(g->getNameOfNode(leftEnd).c_str());	//actual starting position of the node
	int ln = g->getLenOfNode(leftEnd);
	int flag = 0;
	for (int k = 0; k < g->getSizeofGraph(); k++) {
		int tmpSp = atoi(g->getNameOfNode(k).c_str( ));
		int tmpLn = g->getLenOfNode(k);
		if ((sp > tmpSp) && (sp < (tmpSp + tmpLn - 1)) && ((sp + ln - 1)
				> (tmpSp + tmpLn - 1))) {
			//if overlap length is less than windowsize, we consider they're not overlapping.
			//if not, we see them as overlapping,so this is not a real left end.
			if ((tmpSp + tmpLn - sp) >= ((g->ovl).d2.getWindowSize())) {
				out << ret << "Node " << leftEnd << " is not a real left-most node. ";
				ret = out.str();
				out.str("");
				out.clear();
				out << ret << "Overlap length with node " << k << " is " << (tmpSp + tmpLn - sp) << "\n"; //(sp+ln-tmpSp));
				ret = out.str();
				out.str("");
				out.clear();
				flag = 1;
			}
		}
	}

	if (flag == 0) {
		out << ret << "Node " << leftEnd << " is a real left-most node.\n";
		ret = out.str();
	}
	return ret;
}
