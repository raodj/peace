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
	idxOfSNPFile = 0;
}

Reconstruction::Reconstruction(Graph* graph, vector<SixTuple*> align, vector<SixTuple*> leftEnds, InclusionNodes* inc, const std::string& con, const std::string& sing, const std::string& numF, int longestEstLen) {
	g = graph;
	alignArray = align;
	incNodes = inc;
	leftMostNodes = leftEnds;
	consensusFileName = con;
	singletonFileName = sing;
	numOfUsedESTsFileName = numF;
	COMPARISON_LENGTH = longestEstLen;
	idxOfSNPFile = 0;
	idxOfContig = 1; //starting from 1.
	totalNumOfRD = 0;

	usedNodes = vector<int> (g->graphNodes.size(), 0);

	alignment.setScoringSystem(2, -1, -1);
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
    int index = 1;
	vector<string> consensus = reconstruct();
    if (OUTPUT_ACE == 0) {
        outFile1.open(consensusFileName.c_str(),ios::trunc);
        for (size_t i=0; i<consensus.size(); i++) {
            string str = consensus[i];
            size_t found = str.find('\n');

            if (found != string::npos) { //there is "\n" in the sequence
                outFile1 << ">Contig" << (index++) << endl;

                int start = 0;
                while (found != string::npos) {
                    outFile1 << str.substr(start, (int)found-start) << endl;
                    start = found + 1;
                    found = str.find('\n', found+1);
                }
                outFile1 << str.substr(start) << endl;
            } else {
                outFile1 << ">Contig" << (index++) << endl;
                outFile1 << consensus[i] << endl;
            }
        }
        outFile1.close();
    }
    
    if (OUTPUT_ACE == 1) { //output an ACE file with the name "$consensusFileName.ace"
        string tmpFileName = "tmp_" + consensusFileName;
        ifstream ifile(tmpFileName.c_str());
        if (ifile) {
            string aceFileName = consensusFileName; // + ".ace";
            outFile1.open(aceFileName.c_str(),ios::trunc);
            outFile1 << "AS " << index-1 << " " << totalNumOfRD << endl << endl;
            while (ifile.good()) {
                string str;
                getline(ifile, str);
                outFile1 << str << endl;
            }
            outFile1.close();
        }
        ifile.close();
        remove(tmpFileName.c_str());
    }


    /*
     * print singletons and numOfUsedESTs
     */
    ofstream outFile2;
    outFile2.open(singletonFileName.c_str(),ios::trunc);

    int num = 0;
    for (size_t i=0; i<usedNodes.size(); i++) {
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

    vector<string> ret = reconstuctForEachCluster();

    //for debug
    std::stringstream out;
    out << printStr << ret.size() << " consensus from above " << leftMostNodes.size() << " left ends:\n";
    printStr = out.str();
    for (size_t p=0; p<ret.size(); p++) {
        printStr = printStr + ret[p] + "\n";
    }

    //print debug information about the generated consensus.
    cout << "*********************consensus:********************\n";
    cout << printStr << endl;

    return ret;
}

/*
 * Called by reconstruct().
 * Do reconstruction for all the left ends in each cluster. Then return all the consensus sequences.
 */
std::vector<std::string> Reconstruction::reconstuctForEachCluster() {
    vector<vector<int> > dGraph = genDGraph();

    vector<string> retStr;
    vector<int> allNodes(g->graphNodes.size(), 0);

    //assign each node a cluster index and store it into a map
    int clusterIdx = 0;
    map<int, int> curMSTNodes;
    int first = 0;
    while (true) {
        bool exit = 1;
        for (size_t i=first; i<allNodes.size(); i++) {
            if (allNodes[i] == 0) { //find the first node which has not been processed
                first = i;
                allNodes[first] = 1;
                exit = 0;
                break;
            }
        }
        if (exit) break;

        curMSTNodes[first] = clusterIdx;
        int sizeOfCurCluster = 1;

        vector<stack<int> > nodes (2);
        nodes[0].push(first);
        nodes[1].push(-1);
        while (true) { //put all the nodes in the same MST with 'first' into curMSTNodes
            nodes = g->getNodesFromMST(nodes);
            if (nodes[0].size() == 0) break;
            stack<int> tmp = nodes[0];
            while (!tmp.empty()) {
                int tmpIndex = tmp.top();
                tmp.pop();
                curMSTNodes[tmpIndex] = clusterIdx;
                allNodes[tmpIndex] = 1;
                sizeOfCurCluster++;
            }
        }

        if (sizeOfCurCluster <= 1) { //first is a singleton.
            curMSTNodes[first] = -1;
            clusterIdx -= 1;
        }

        clusterIdx++;
        first++;
        if (sizeOfCurCluster > 1)
            cout << "Cluster " << clusterIdx-1 << " , size is " << sizeOfCurCluster << endl;
    }

    //put those left ends in one cluster into a vector
    vector<int>  clusterLeftEnds[clusterIdx];
    int sizeOfs = leftMostNodes.size();
    for (int i=0; i<sizeOfs; i++) { //start for
        int idx = leftMostNodes[i]->curNode;
        map<int, int>::iterator ite = curMSTNodes.find(idx);
        if ((ite != curMSTNodes.end()) && (ite->second != -1)) { //idx is found and is not a singleton
            clusterLeftEnds[ite->second].push_back(idx);
        }
    } //end for

    //call processLeftEnds to process these left ends cluster by cluster
    for (int i=0; i<clusterIdx; i++) {
        cout << "Cluster " << i << ", number of left ends in the cluster is " << clusterLeftEnds[i].size() << endl;
        vector<string> ret = processLeftEnds(clusterLeftEnds[i], dGraph); //reconstruct from left ends in the cluster
        for (size_t k=0; k<ret.size(); k++)
            retStr.push_back(ret[k]);
    }

    return retStr;
}


std::vector<std::vector<int> > Reconstruction::genDGraph() {
    /*
     * Calculate the length of dGraph.
     */
    int len = 0;
    for (size_t i=0; i<alignArray.size(); i++) {
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
    for (size_t i=0; i<alignArray.size(); i++) {
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

    // Make a directed graph graphForPrim.
    graphForPrim = DefGraph(g->graphNodes.size());
    for (size_t j=0; j<dGraph.size(); j++) {
        if (dGraph[j][3] != 0) {	//there is an edge between the nodes
            graphForPrim[dGraph[j][0]].push_back(Edge(dGraph[j][1], dGraph[j][2]));
        }
    }


    return dGraph;
}



/*
 * This method is designed to group left ends into independent sets. All the dependent left ends (include one by one) are put
 * into one group.
 * In each group, we put all the nodes which can be connected by any left end in the group,
 * and combine them together to form a consensus by calling the function "processLeftEndsWithInclusion".

 * @return all the consensus sequences.
 */
vector<string>Reconstruction::processLeftEnds(vector<int>& curLeftEnds, vector<vector<int> >& dGraph) {
    vector<string> allOutputContigs;	//store all the generated sequences
    int sizeOfs = curLeftEnds.size();
    if (sizeOfs == 0) {
        cout << "zero left ends in the cluster!" << endl;
        allOutputContigs.push_back("");
        return allOutputContigs;
    }
    vector<LeftEnd> resultArray; //store the starting positions of ests
    for (int i=0; i<sizeOfs; i++) { //start for
        int idx = curLeftEnds[i];
        resultArray.push_back(LeftEnd(idx, g->getSeqOfNode(idx)));
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
        for (size_t i=1; i<allLeftEnds.size(); i++) {
            vector<int> ovlDis = g->calDist.searchDistance(allLeftEnds[i].index, s1Idx);

            if (ovlDis[1] == INT_MIN) { //if resultArray[i].firstEst is included in s1
                includedEnds.push_back(allLeftEnds[i]);
                incNodes->addNode2(allLeftEnds[i].index, s1Idx);
            } else if ((ovlDis[1] == INT_MAX) && (ovlDis[0] == INT_MAX)) { //have not been calculated
                bool b;
                vector<int> qual1 = g->graphNodes[allLeftEnds[i].index].getQualScores();
                vector<int> qual2 = g->graphNodes[s1Idx].getQualScores();
                if (i == 1) {
                    if (USE_QUALITY_FILE == 1) { //use quality file
                        b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1, true, qual1, qual2); //not keep hash table of s1 in D2.cpp
                    } else {
                        b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1, true); //not keep hash table of s1 in D2.cpp
                    }
                } else {
                    if (USE_QUALITY_FILE == 1) { //use quality file
                        b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1, false, qual1, qual2); //not keep hash table of s1 in D2.cpp
                    } else {
                        b = (g->ovl).checkInclusion(allLeftEnds[i].seq, s1, false);//keep hash table of s1 in D2.cpp
                    }
                }

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
 * We put all the nodes together which can be connected by any left end in the group,
 * and reconstruct a consensus from them.
 */
string Reconstruction::processLeftEndsWithInclusion(vector<LeftEnd>& includeStrs, vector<vector<int> >& dGraph) {
    cout << includeStrs.size() << endl;
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
            AlignResult strs;
            if (USE_QUALITY_FILE == 1) { //use quality file
                vector<int> qual1 = g->getQualScoresOfNode(includeStrs[idxMaxLen].index);
                vector<int> qual2 = g->getQualScoresOfNode(includeStrs[i].index);
                strs = alignment.getLocalAlignment(s1, s2, qual1, qual2);
            } else {
                strs = alignment.getLocalAlignment(s1, s2);
            }
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
 *  Return consensus from this set of left ends.
 *  input: leftEnds[][2] : [][0]-index, [][1]-position of the left end; size-the size of the array
 *
 * For each left-end node, starting from it to calculate positions for each node.
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
        DefGraph primMST = constructMinTree(leftEnd); //the first param is the total number of ESTs.

        //get starting positions for the nodes in primMST
        getStartPos(leftEnd, primMST, dGraph);

        //printStr = printStr + printLeftEndInfo(leftEnd); //get information of this left end which is used to start reconstruction.
    }
    stringstream tstr;
    tstr << leftEnds[0][0];
    printStr = printStr + "Start from node " + tstr.str() + " to do reconstruction" + "\n\n";

    vector<StartPos> tmpArray;
    for (size_t j=0; j<sPos.size(); j++) {
        ite = ends.find(j);
        if ((ite != ends.end()) || (sPos[j] != 0)) { // is left end or a node which has got its position
            tmpArray.push_back(StartPos(sPos[j], j));
        }
    }


    vector<string> tStr;
    if (USE_QUALITY_FILE == 1) { //use quality file
        tStr = reconstructSeqWithQual(tmpArray);
    } else {
        tStr = reconstructSeq(tmpArray);
    }
    if (tStr.size() != 0) {
        printStr = printStr + tStr[1] + " nodes are used to reconstruct the sequence.\n";
        printStr = printStr + tStr[0] + "\n\n";
        ret = tStr[0];
    }

    return ret;
}

/*
 * reconstruct a sequence which starts from a left end.
 * return: ret[0]-the consensus, ret[1]-the total number of nodes used for reconstruction.
 */
vector<string> Reconstruction::reconstructSeq(vector<StartPos>& a) {
    std::stringstream out;
    vector<string> ret(2);
    if (a.size() == 0) {
        ret[1] = "0";
        return ret;
    }

    sort(a.begin(), a.end());
    vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction and set usedNodes[related nodes]=1.
    out << addedNodes.size();
    ret[1] = out.str();
    if (addedNodes.size() == 1){
        //size=1, this node should be a singleton if it is not used in other place. it shouldn't appear as a consensus sequence.
        //size=1, then usedNodes[this node] will not be set to 1 in addInclusionNodes method.
        return ret;
    }
    vector<StartPos> resultArray(addedNodes.size());
    for (size_t i=0; i<addedNodes.size(); i++) {
        UsedNode tmpNode = addedNodes[i];
        resultArray[i] = StartPos(tmpNode.pos, tmpNode.index);
    }
    sort(resultArray.begin(), resultArray.end());

    //Begin: For ACE output format
    string COStr = "";
    string BQStr = "BQ\n";
    string AFStr = "";
    string BSStr = "";
    string RDStr = "";
    int curBS = 0;
    int numOfBS = 0;
    std::stringstream ss;
    //End: For ACE output format

    int comparisonLen = COMPARISON_LENGTH;
    vector<SingleBase*> bases;
    string tConsensus = g->getSeqOfNode(resultArray[0].index);
    string Name0 = g->getNameOfNode(resultArray[0].index);
    string curSeq = "";
    int len = resultArray.size() - 1;
    for (int i=1; i<=len; i++) {
        curSeq = g->getSeqOfNode(resultArray[i].index);
        string curName = g->getNameOfNode(resultArray[i].index);
        string tmpConsensus = tConsensus;
        if ((int) tmpConsensus.length() > comparisonLen) {
            tmpConsensus = tConsensus.substr(tConsensus.size()-comparisonLen+1);
        }

        AlignResult strs;
        if (USE_BOUNDED_SW == 0) { //use ordinary version
            strs = alignment.getLocalAlignment(tmpConsensus, curSeq);
        } else  { //use bounded version
            strs = alignment.getBoundedLocalAlignment(tmpConsensus, curSeq);
        }
        if ((strs.str1.size()==0) || (strs.str2.size()==0)) continue;
        tConsensus = replace(tConsensus, replace(strs.str1, "-", ""), strs.str1);
        int offset = (int)tConsensus.find(strs.str1);
        if (offset < 0) continue; //not found
        string tSeq = replace(curSeq, replace(strs.str2, "-", ""), strs.str2);
        int tmpOff = (int)tSeq.find(strs.str2);
        if (tmpOff < 0) continue; //not found
        string firstPartCur = tSeq.substr(0, tmpOff);
        curSeq = tSeq.substr(tmpOff);
        if (i == 1) {
            //Begin: For ACE output format
            if (OUTPUT_ACE == 1) {
                if (offset > tmpOff) { //use the first read as the first part in the consensus
                    AFStr += "AF " + Name0 + " " + g->getDirectionOfNode(resultArray[0].index) + " 1" + "\n";
                    ss << "AF " << curName << " " << g->getDirectionOfNode(resultArray[i].index) <<" "<< offset-tmpOff+1 << "\n";
                    AFStr += ss.str();
                    ss.str("");
                    ss.clear();

                    ss << "BS 1 " << offset-tmpOff <<  " " << Name0 << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();
                    curBS = offset+strs.str1.size();
                    ss << "BS " << offset-tmpOff+1 << " " << curBS << " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();

                    numOfBS = 2;
                } else {
                    AFStr += "AF " + curName + " " + g->getDirectionOfNode(resultArray[0].index) + " 1" + "\n";
                    ss << "AF " << Name0 << " " << g->getDirectionOfNode(resultArray[i].index) <<" "<< tmpOff-offset+1 << "\n";
                    AFStr += ss.str();
                    ss.str("");
                    ss.clear();

                    curBS = tmpOff+strs.str2.size();
                    ss << "BS 1 " << curBS <<  " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();

                    numOfBS = 1;
                }
                ss << tConsensus.size();
                RDStr += "RD " + Name0 + " " + ss.str() + " 0 0\n";
                RDStr += replace(tConsensus, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();

                ss << tSeq.size();
                RDStr += "RD " + curName + " " + ss.str() + " 0 0\n";
                RDStr += replace(tSeq, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();
            }
            //End: For ACE output format

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
            int end = offset+len2;

            //Begin: For ACE output format
            if (OUTPUT_ACE == 1) {
                ss << offset-tmpOff+1;
                AFStr += "AF " + curName + " " + g->getDirectionOfNode(resultArray[i].index) + " " + ss.str() + "\n";
                ss.str("");
                ss.clear();

                int tmpBS = offset + strs.str2.size();
                if (tmpBS > curBS) { //write a new BS item only if tmpBS is bigger than curBS.
                    ss << "BS " << curBS+1 << " " << tmpBS << " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();
                    curBS = tmpBS;

                    numOfBS++;
                }

                ss << tSeq.size();
                RDStr += "RD " + curName + " " + ss.str() + " 0 0\n";
                RDStr += replace(tSeq, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();
            }
            //End: For ACE output format

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

    int totalLen = bases.size();
    /*
      for (int i=bases.size()-1; i>=0; i--) {
      if (bases[i]->getTotalNumOfBase() > 2) {
      break;
      }
      totalLen--;
      }
      tConsensus = tConsensus.substr(0, totalLen);
    */

    ret[0] = replace(tConsensus, "P", "");

    if (OUTPUT_SNP == 1) { //write output files for SNP analysis
        writeToSNPFile(bases, totalLen);
    }

    //Begin: For ACE output format
    if (OUTPUT_ACE == 1) {
        ss << "CO Contig" << idxOfContig << " " << tConsensus.size() << " " << ret[1] << " " << numOfBS << " U\n";
        COStr += ss.str();
        ss.str("");
        ss.clear();
        COStr += replace(tConsensus, "P", "*") + "\n";
        idxOfContig++;

        for (size_t i=0; i<bases.size(); i++) {
            int phredScore = bases[i]->getQualScore();
            if (phredScore != -1) {
                ss << phredScore;
                BQStr += " " + ss.str();
                ss.str("");
                ss.clear();
            }
        }
        BQStr += "\n";

        //write to the temporary ACE file
        string tmpFileName = "tmp_" + consensusFileName;
        ofstream outFile;
        outFile.open(tmpFileName.c_str(),ios::app);
        outFile << COStr << endl;
        outFile << BQStr << endl;
        outFile << AFStr << endl;
        outFile << BSStr << endl;
        outFile << RDStr << endl;
        outFile.close();
    }
    //End: For ACE output format

    //release bases
    for (size_t i=0; i<bases.size(); i++) {
        delete bases[i];
    }

    return ret;
}

/*
 * reconstruct a sequence which starts from a left end using quality scores.
 * return: ret[0]-the consensus, ret[1]-the total number of nodes used for reconstruction.
 */
vector<string> Reconstruction::reconstructSeqWithQual(vector<StartPos>& a) {
    std::stringstream out;
    vector<string> ret(2);
    if (a.size() == 0) {
        ret[1] = "0";
        return ret;
    }

    sort(a.begin(), a.end());
    vector<UsedNode> addedNodes = addInclusionNodes(a);  //add all those related inclusion nodes into it for reconstruction and set usedNodes[related nodes]=1.
    out << addedNodes.size();
    ret[1] = out.str();
    if (addedNodes.size() == 1){
        //size=1, this node should be a singleton if it is not used in other place. it shouldn't appear as a consensus sequence.
        //size=1, then usedNodes[this node] will not be set to 1 in addInclusionNodes method.
        return ret;
    }
    vector<StartPos> resultArray(addedNodes.size());
    for (size_t i=0; i<addedNodes.size(); i++) {
        UsedNode tmpNode = addedNodes[i];
        resultArray[i] = StartPos(tmpNode.pos, tmpNode.index);
    }
    sort(resultArray.begin(), resultArray.end());
    //print resultArray for debug
    /*
      for (int i=0; i<resultArray.size(); i++) {
      cout << resultArray[i].index << endl;
      }*/

    //Begin: For ACE output format
    string COStr = "";
    string BQStr = "BQ\n";
    string AFStr = "";
    string BSStr = "";
    string RDStr = "";
    int curBS = 0;
    int numOfBS = 0;
    std::stringstream ss;
    //End: For ACE output format

    int comparisonLen = COMPARISON_LENGTH;
    vector<SingleBase*> bases;
    string tConsensus = g->getSeqOfNode(resultArray[0].index);
    vector<int> firstSeqQualScores = g->getQualScoresOfNode(resultArray[0].index);
    string Name0 = g->getNameOfNode(resultArray[0].index);
    string curSeq = "";
    int len = resultArray.size() - 1;
    for (int i=1; i<=len; i++) {
        curSeq = g->getSeqOfNode(resultArray[i].index);
        string curName = g->getNameOfNode(resultArray[i].index);
        vector<int> curSeqQualScores = g->getQualScoresOfNode(resultArray[i].index);
        string tmpConsensus = tConsensus;

        vector<int> qualOfTmpConsensus;
        if ((int) tmpConsensus.length() > comparisonLen) {
            int tmpConStart = tConsensus.size()-comparisonLen+1;
            tmpConsensus = tConsensus.substr(tmpConStart);
            if (i == 1) {
                for (size_t tmpCon=tmpConStart; tmpCon < tConsensus.size(); tmpCon++) {
                    qualOfTmpConsensus.push_back(firstSeqQualScores[tmpCon]);
                }
            } else {
                for (size_t tmpCon=tmpConStart; tmpCon < tConsensus.size(); tmpCon++) {
                    qualOfTmpConsensus.push_back(bases[tmpCon]->curqualVal);
                }
            }
        } else if (i == 1) {
            qualOfTmpConsensus = firstSeqQualScores;
        } else {
            for (size_t tmpCon=0; tmpCon < tConsensus.size(); tmpCon++) {
                qualOfTmpConsensus.push_back(bases[tmpCon]->curqualVal);
            }
        }

        AlignResult strs;
        if (USE_BOUNDED_SW == 0) { //use ordinary version
            strs = alignment.getLocalAlignment(tmpConsensus, curSeq, qualOfTmpConsensus, curSeqQualScores);
        } else  { //use bounded version
            strs = alignment.getBoundedLocalAlignment(tmpConsensus, curSeq);
        }
        if ((strs.str1.size()==0) || (strs.str2.size()==0)) continue;
        tConsensus = replace(tConsensus, replace(strs.str1, "-", ""), strs.str1);
        int offset = (int)tConsensus.find(strs.str1);
        if (offset < 0) continue; //not found

        string tSeq = replace(curSeq, replace(strs.str2, "-", ""), strs.str2);
        int tmpOff = (int)tSeq.find(strs.str2);
        if (tmpOff < 0) continue; //not found
        int posInEst = tmpOff;
        string firstPartCur = tSeq.substr(0, tmpOff);
        curSeq = tSeq.substr(tmpOff);
        if (i == 1) {
            //Begin: For ACE output format
            if (OUTPUT_ACE == 1) {
                if (offset > tmpOff) { //use the first read as the first part in the consensus
                    AFStr += "AF " + Name0 + " " + g->getDirectionOfNode(resultArray[0].index) + " 1" + "\n";
                    ss << "AF " << curName << " " << g->getDirectionOfNode(resultArray[i].index) <<" "<< offset-tmpOff+1 << "\n";
                    AFStr += ss.str();
                    ss.str("");
                    ss.clear();

                    ss << "BS 1 " << offset-tmpOff <<  " " << Name0 << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();
                    curBS = offset+strs.str1.size();
                    ss << "BS " << offset-tmpOff+1 << " " << curBS << " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();

                    numOfBS = 2;
                } else {
                    AFStr += "AF " + curName + " " + g->getDirectionOfNode(resultArray[0].index) + " 1" + "\n";
                    ss << "AF " << Name0 << " " << g->getDirectionOfNode(resultArray[i].index) <<" "<< tmpOff-offset+1 << "\n";
                    AFStr += ss.str();
                    ss.str("");
                    ss.clear();

                    curBS = tmpOff+strs.str2.size();
                    ss << "BS 1 " << curBS <<  " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();

                    numOfBS = 1;
                }
                ss << tConsensus.size();
                RDStr += "RD " + Name0 + " " + ss.str() + " 0 0\n";
                RDStr += replace(tConsensus, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();

                ss << tSeq.size();
                RDStr += "RD " + curName + " " + ss.str() + " 0 0\n";
                RDStr += replace(tSeq, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();
            }
            //End: For ACE output format

            int len1 = tConsensus.size();
            int len2 = curSeq.size();
            int end = max(offset+len2, len1);
            int curScorePosinS1 = 0;
            int curScorePosinS2 = 0;
            for (int j=0; j<end; j++) {
                int qualScore1 = 0;
                int qualScore2 = 0;

                if ((j < len1) && (j-offset >= 0) && (j-offset < len2)) { //overlap part
                    char c1 = tConsensus[j];
                    char c2 = curSeq[j-offset];
                    if (c1 != '-') {
                        qualScore1 = firstSeqQualScores[curScorePosinS1];
                        curScorePosinS1++;
                    } else {
                        qualScore1 = curSeqQualScores[posInEst+curScorePosinS2]; //equal to qualScore2
                    }
                    if (c2 != '-') {
                        qualScore2 = curSeqQualScores[posInEst+curScorePosinS2];
                        curScorePosinS2++;
                    } else {
                        qualScore2 = qualScore1; //equal to qualScore1
                    }
                    bases.push_back(new SingleBase(c1, c2, qualScore1, qualScore2));
                } else if ((j-offset < 0) || (j-offset >= len2)) {
                    char c1 = tConsensus[j];
                    if (c1 != '-') {
                        qualScore1 = firstSeqQualScores[curScorePosinS1];
                        curScorePosinS1++;
                    }
                    bases.push_back(new SingleBase(c1, qualScore1));
                } else if (j >= len1) {
                    char c2 = curSeq[j-offset];
                    if (c2 != '-') {
                        qualScore2 = curSeqQualScores[posInEst+curScorePosinS2];
                        curScorePosinS2++;
                    }
                    bases.push_back(new SingleBase(c2, qualScore2));
                }
            }
        } else {
            int len1 = tConsensus.size();
            int len2 = curSeq.size();
            int end = offset+len2;
            int curScorePos = 0;

            //Begin: For ACE output format
            if (OUTPUT_ACE == 1) {
                ss << offset-tmpOff+1;
                AFStr += "AF " + curName + " " + g->getDirectionOfNode(resultArray[i].index) + " " + ss.str() + "\n";
                ss.str("");
                ss.clear();

                int tmpBS = offset + strs.str2.size();
                if (tmpBS > curBS) { //write a new BS item only if tmpBS is bigger than curBS.
                    ss << "BS " << curBS+1 << " " << tmpBS << " " << curName << "\n";
                    BSStr += ss.str();
                    ss.str("");
                    ss.clear();
                    curBS = tmpBS;

                    numOfBS++;
                }

                ss << tSeq.size();
                RDStr += "RD " + curName + " " + ss.str() + " 0 0\n";
                RDStr += replace(tSeq, "-", "*") + "\n\n";
                RDStr += "QA 1 " + ss.str() + " 1 " + ss.str() + "\n\n";
                totalNumOfRD++;
                ss.str("");
                ss.clear();
            }
            //End: For ACE output format

            for (int j=offset; j<end; j++) {
                char c2 = curSeq[j-offset];
                int qualScore = 0;
                if (c2 != '-') {
                    qualScore = curSeqQualScores[posInEst+curScorePos];
                    curScorePos++;
                } else if ((int) bases.size() > j){
                    int curBaseNum = bases[j]->curBaseNum;
                    if (curBaseNum > 0) {
                        qualScore = (bases[j]->curqualVal)/curBaseNum;
                    }
                }

                if ((j < len1) && (j-offset < len2)) { //overlap part
                    char c1 = tConsensus[j];

                    if (c1 != '-') {
                        bases[j]->addOneBase(c2, qualScore);
                    } else {
                        bases.insert(bases.begin()+j, new SingleBase(c1, c2, qualScore, qualScore));
                    }
                } else if (j >= len1) {
                    bases.push_back(new SingleBase(c2, qualScore));
                }
            }

            //put first part of curSeq into bases
            for (int k=0; k<tmpOff; k++) {
                if (offset-tmpOff+k >= 0)
                    bases[offset-tmpOff+k]->addOneBase(firstPartCur[k], curSeqQualScores[k]);
            }
        }

        tConsensus = getCurConsensus(bases);
    }

    int totalLen = bases.size();
    /*
      for (int i=bases.size()-1; i>=0; i--) {
      if (bases[i]->getTotalNumOfBase() > 2) {
      break;
      }
      totalLen--;
      }
      tConsensus = tConsensus.substr(0, totalLen);
    */

    ret[0] = replace(tConsensus, "P", "");

    if (OUTPUT_SNP == 1) { //write output files for SNP analysis
        writeToSNPFile(bases, totalLen);
    }

    //Begin: For ACE output format
    if (OUTPUT_ACE == 1) {
        ss << "CO Contig" << idxOfContig << " " << tConsensus.size() << " " << ret[1] << " " << numOfBS << " U\n";
        COStr += ss.str();
        ss.str("");
        ss.clear();
        COStr += replace(tConsensus, "P", "*") + "\n";
        idxOfContig++;

        for (size_t i=0; i<bases.size(); i++) {
            int phredScore = bases[i]->getQualScore();
            if (phredScore != -1) {
                ss << phredScore;
                BQStr += " " + ss.str();
                ss.str("");
                ss.clear();
            }
        }
        BQStr += "\n";

        //write to the temporary ACE file
        string tmpFileName = "tmp_" + consensusFileName;
        ofstream outFile;
        outFile.open(tmpFileName.c_str(),ios::app);
        outFile << COStr << endl;
        outFile << BQStr << endl;
        outFile << AFStr << endl;
        outFile << BSStr << endl;
        outFile << RDStr << endl;
        outFile.close();
    }
    //End: For ACE output format

    //release bases
    for (size_t i=0; i<bases.size(); i++) {
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

void Reconstruction::writeToSNPFile(vector<SingleBase*> bases, int totalLen) {
    if (totalLen <= 0) return;
    stringstream ss;
    ss << idxOfSNPFile++;
    string outFileName = "snpAnalysis." + ss.str();
    ofstream outFile1;
    outFile1.open(outFileName.c_str(),ios::trunc);

    outFile1 << "A\tG\tC\tT" << endl;
    for (int i=0; i<totalLen; i++) {
        outFile1 << bases[i]->numA << "\t" << bases[i]->numG << "\t" << bases[i]->numC << "\t" << bases[i]->numT << endl;
    }
    outFile1.close();
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

    if (old.size() == 0) return str; //old is ""
    int newSize = newstr.size();
    size_t found = str.find(old);
    while (found != string::npos) { //there is "old" in the sequence
        str.replace(found, old.size(), newstr);
        if ((found+newSize)<str.size())
            found = str.find(old, found+newSize);
        else
            found = string::npos;
    }
    return str;

    /*
      size_t found = str.find(old);
      while (found != string::npos) { //there is "old" in the sequence
      str.replace(found, old.size(), newstr);
      found = str.find(old, found+old.size()-1);
      }
      return str;
    */

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
        for (size_t j=0; j<chdIdx.size(); j++) {
            tmpList.insert(pair<int, int> (chdIdx[j], pos+1));
            usedNodes[chdIdx[j]] = 1; //mark this node being used
        }
    }
    for (map<int, int>::iterator it = tmpList.begin(); it != tmpList.end(); it++) {
        retList.push_back(UsedNode(it->first, it->second));
    }
    sort(retList.begin(), retList.end());

    if (retList.size()==1) //if only one node in retList after adding inclusion nodes, this node should be a singleton.
        usedNodes[retList[0].index] = 0;
    return retList;
}

/*
 * Construct a directed Miminum spanning tree.
 *
 *  @param nOfNodes number of nodes
 *  @param g a directed graph, the second dimension has three elements:
 *  	index of starting node, index of ending node, weight between them.
 */
DefGraph Reconstruction::constructMinTree(int source) {
    DefGraph mst = Prim(graphForPrim, source, true);
    return mst;
}

/*
 * Calculate starting positions for each node.
 */
void Reconstruction::getStartPos(int parentNode, DefGraph& tree, vector<vector<int> >& d) {
    for (EdgeIterator it=tree[parentNode].begin(); it!=tree[parentNode].end(); it++) {
        int index = (*it).node;

        int overlapLen = 0;
        for (size_t i=0; i<d.size(); i++) {
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
