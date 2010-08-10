#include "ProcessRC.h"
#include <map>
#include <cstdlib>
using namespace std;


ProcessRC::ProcessRC(const std::string& estF, const std::string& mstF, const std::string& newEstF) {
    estFileName = estF;
    mstFileName = mstF;
    newEstFileName = newEstF;
    treeThreshold = 1;

	string aceInfoFile = estFileName + ".aceinfo";
	aceOutFile.open(aceInfoFile.c_str(), ios::trunc);
}

void ProcessRC::operateRC() {
	cout << "start readEstFile" << endl;
	readEstFile(estFileName);
	cout << "start readMST" << endl;
	readMST(mstFileName);
	cout << "start ProcessEachCluster" << endl;
	ProcessEachCluster();
	cout << "start writeEstFile" << endl;
	writeEstFile();
	aceOutFile.close();
}

/*
 * Read ests from the input file;
 * Generate a Graph object, all the ests are considered to be one node in the graph;
 * No edge in the graph. Edges will be added in "createAlignArray" function.
 *
 * FASTA format:
 * A sequence in FASTA format begins with a single-line description, followed by lines of
 * sequence data. The description line is distinguished from the sequence data by a greater-than
 * (">") symbol in the first column. The word following the ">" symbol is the identifier of the
 * sequence, and the rest of the line is the description (both are optional). There should be no
 * space between the ">" and the first letter of the identifier. It is recommended that all lines
 * of text be shorter than 80 characters. The sequence ends if another line starting with a ">"
 * appears; this indicates the start of another sequence.
 *
 * EST file comment format:
 * >g001_001169_001679: first one is number of gene, second is index of tarting position, third is
 * 						index of ending position. Index starts from 0.
 */
void ProcessRC::readEstFile(const string& inFileName) {
	vector<string> ests; //store all the ests.

	ifstream in;
	in.open(inFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);
		while (!in.eof()) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str.size()==0) getline(in, str);
			// first line is comment line which begins from '>'
			if (str[0] == '>') {	//comment line begins from '>'
				ests.push_back("0");
				ests.push_back(str.substr(1)); //comment

				//get est in the next lines
				getline(in, str);
				string estStr;
				while (str.size() != 0) {
					str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
					if (str.size() != 0)	{
						if (str[0] != '>') {
							estStr.append(str);
						} else  {
							ests.push_back(estStr);
							break;
						}
					}
					getline(in, str);
				}
				if (str.size() == 0) {
					ests.push_back(estStr);
				}
			}
		}
	}
	in.close();

	//generate a graph from the input file
	int i=0;
	cout << "number of ests: " << ests.size()/3 << endl;
	while (i<ests.size()) {
		//the sequence is upper-case
		graphNodes.push_back(Node(ests[i], ests[i+1], toUpperCase(ests[i+2])));
		i = i+3;
	}
}

/*
 * read a minimum spanning tree from the input MST file
 */
void ProcessRC::readMST(const string& inFileName) {
	int nOfNodes = graphNodes.size();
	if (nOfNodes == 0) {
		cout << "zero nodes in the peace-generated input MST file!" << endl;
		exit(1);
	}

	//vector<vector<int> > nodes(nOfNodes-1, vector<int>(3));	//store edges in MST, there are n-1 edges, n is number of nodes.
	vector<vector<float> > nodes;

	//read mst from the input file
	ifstream in;
	in.open(inFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);

		while (str.size() != 0) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str[0] != '#') {	//comment line begins from '#'
				vector<string> paras = split(str, ',');
				if (atoi(paras[0].c_str()) != -1) {	//-1 means root of MST
					vector<float> tmp(4);
					tmp[0] = atof(paras[0].c_str());
					tmp[1] = atof(paras[1].c_str());
					tmp[2] = atof(paras[2].c_str());
					tmp[3] = atof(paras[4].c_str());
					nodes.push_back(tmp);
				}
			}
			getline(in, str);
		}
	}
	in.close();

	// Make a undirected MST.
	mst = DefGraph(nOfNodes);
	for (int j=0; j<nodes.size(); j++) {
		addEdge(mst, int(nodes[j][0]), int(nodes[j][1]), nodes[j][2], int(nodes[j][3]), false);
	}
}

void ProcessRC::ProcessEachCluster() {
	vector<int> allNodes(graphNodes.size(), 0);
	int maxClusterSize = 0;
	int totalSize = 0;

	while (true) {
		int first = -1;
		bool exit = 1;
		for (int i=0; i<allNodes.size(); i++) {
			if (allNodes[i] == 0) { //find the first node which has not been processed
				first = i;
				allNodes[first] = 1;
				exit = 0;
				break;
			}
		}
		if (exit) break;

		vector<stack<int> > nodes (3);
		nodes[0].push(first);
		nodes[1].push(-1);
		nodes[2].push(1);

		//cout << "\t\t"<<graphNodes[first].comment << endl;
		int size = 1;
		while (true) { //put all the nodes in the same MST with 'first' into curMSTNodes
			nodes = getNodesFromMST(nodes);
			size += nodes[0].size();
			if (nodes[0].size() == 0) break;
			stack<int> tmp = nodes[0];
			while (!tmp.empty()) {
				int tmpIndex = tmp.top();
				tmp.pop();
				allNodes[tmpIndex] = 1;
				//cout << "\t\t"<<graphNodes[tmpIndex].comment << endl;
			}
		}
		cout << "\tsize of the cluster is " << size << endl;
		cout << "\t\tfirst it the cluster is " << first << endl;
		maxClusterSize = (maxClusterSize>size)? maxClusterSize : size;
		totalSize += size;
	}
	cout << "The maximal cluster size is " << maxClusterSize << endl;
	cout << "Total number of nodes is " << graphNodes.size() << endl;
	cout << "Total number of nodes in all the clusters is " << totalSize << endl;
}

vector<stack<int> > ProcessRC::getNodesFromMST(vector<stack<int> > nodes) {
	vector<stack<int> > ret (3);
	while (!(nodes[0].empty())) {
		int curIndex = nodes[0].top();
		aceOutFile << curIndex << " 1\n";
		nodes[0].pop();
		int parentIndex = nodes[1].top();
		nodes[1].pop();
		int reverseFlag = nodes[2].top();
		nodes[2].pop();

		for (EdgeIterator it=mst[curIndex].begin(); it!=mst[curIndex].end(); it++) {
			int index2 = (*it).node;
			float weight = (*it).weight;
			if (weight <= treeThreshold) {
				if (index2 != parentIndex) {
					ret[0].push(index2);
					ret[1].push(curIndex);
					reverseFlag = reverseFlag * ((*it).reverse);
					aceOutFile << index2 << " " << reverseFlag << "\n";
					ret[2].push(reverseFlag);
					if (reverseFlag < 0) { //reverse the corresponding sequence
						string str = graphNodes[index2].sequence;
						//cout << "ori str: " << str << endl;
						string newStr = "";
						for (int i=0; i<str.size(); i++) {
							switch (str.at(i)) {
								case 'A':
									newStr += "T";
									break;
								case 'T':
									newStr += "A";
									break;
								case 'C':
									newStr += "G";
									break;
								case 'G':
									newStr += "C";
									break;
							}
						}
						reverse(newStr.begin(), newStr.end());
						//cout << "new str: " << str << endl << endl;
						graphNodes[index2].sequence = newStr;
					}
				}
			}
		}
	}
	return ret;
}

vector<string> ProcessRC::split(const string& str, char delimit) {
	vector<string> retStr;
	size_t found = str.find(delimit);
	if (found != string::npos) {
		int start = 0;
		while (found != string::npos) {
			retStr.push_back(str.substr(start, (int) found - start));
			start = found + 1;
			found = str.find(delimit, found + 1);
		}
		retStr.push_back(str.substr(start));
	}
	return retStr;
}

string ProcessRC::toUpperCase(const string& str) {
	string retStr;
	for (int i=0; i<str.size(); i++) {
		retStr.push_back(toupper(str[i]));
	}
	return retStr;
}

void ProcessRC::writeEstFile() {
	ofstream outFile;
	outFile.open(newEstFileName.c_str(), ios::trunc);

	for (int i=0; i<graphNodes.size(); i++) {
		outFile << ">" << graphNodes[i].comment << endl;
		outFile << graphNodes[i].sequence << endl;
	}
	outFile.close();
}

void ProcessRC::separateCluster(int idx, const std::string& newEstFile) {
	int first = idx;
	vector<int> nodesInCluster;
	nodesInCluster.push_back(first);

	vector<stack<int> > nodes (3);
	nodes[0].push(first);
	nodes[1].push(-1);
	nodes[2].push(1);
	int size = 1;
	while (true) { //put all the nodes in the same MST with 'first' into curMSTNodes
		nodes = getNodesFromMST(nodes);
		size += nodes[0].size();
		if (nodes[0].size() == 0) break;
		stack<int> tmp = nodes[0];
		while (!tmp.empty()) {
			int tmpIndex = tmp.top();
			tmp.pop();
			nodesInCluster.push_back(tmpIndex);
		}
	}
	cout << "\tsize of the cluster is " << size << endl;

	ofstream outFile;
	outFile.open(newEstFile.c_str(), ios::trunc);

	for (int i=0; i<nodesInCluster.size(); i++) {
		outFile << ">" << graphNodes[nodesInCluster[i]].comment << endl;
		outFile << graphNodes[nodesInCluster[i]].sequence << endl;
	}
	outFile.close();

}
