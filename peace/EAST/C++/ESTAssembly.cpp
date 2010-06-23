#include "ESTAssembly.h"
#include <cstdlib>
using namespace std;


ESTAssembly::ESTAssembly(std::string& estF, std::string& mstF) {
	incNodes = new InclusionNodes();
	g = new Graph(incNodes);
    gen = NULL;
    rec = NULL;
    estFileName = estF;
    mstFileName = mstF;
}

ESTAssembly::ESTAssembly(std::string& estF, std::string& mstF, std::string& qualF) {
	incNodes = new InclusionNodes();
	g = new Graph(incNodes);
    gen = NULL;
    rec = NULL;
    estFileName = estF;
    mstFileName = mstF;
    qualFileName = qualF;
}

void ESTAssembly::assemble(std::string& con, std::string& sing, std::string& numF) {
	int longestEstLen = 0;
	if (USE_QUALITY_FILE == 1) { //use quality file
		longestEstLen = readEstQualFile();
	} else {
		longestEstLen = readEstFile();
	}
	readMST();
	gen = new SixTuplesGeneration(g, incNodes);
	rec = new Reconstruction(g, gen->getAlignArray(), gen->getLeftEnds(), incNodes, con, sing, numF, longestEstLen);
	rec->getConsensus();
}

ESTAssembly::~ESTAssembly() {
	delete incNodes;
	delete g;
	delete gen;
	delete rec;
}
/*
 * Read ests from the input file;
 * Generate a Graph object, all the ests are considered to be one node in the graph;
 * No edge in the graph. Edges will be added in "createAlignArray" function.
 *
 * Return the length of the longest EST which will be used as COMPARISON_LENGTH in reconstruction.
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
int ESTAssembly::readEstFile() {
	vector<string> ests; //store all the ests.
	int longestLen = 0;

	ifstream in;
	in.open(estFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);
		while (!in.eof()) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str.size()==0) getline(in, str);
			// first line is comment line which begins from '>'
			if (str[0] == '>') {	//comment line begins from '>'
				//vector<string> paras = split(str, '_');
				//ests.push_back(paras[1]);
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
							int theSize = estStr.size();
							if (theSize > longestLen)
								longestLen = theSize;
							break;
						}
					}
					getline(in, str);
				}
				if (str.size() == 0) {
					int theSize = estStr.size();
					if (theSize > longestLen)
						longestLen = theSize;
					ests.push_back(estStr);
				}
			}
		}
	}
	in.close();

	//generate a graph from the input file
	int i=0;
	while (i<ests.size()) {
		//the sequence is upper-case
		g->addNode(Node(ests[i], ests[i+1], toUpperCase(ests[i+2])));
		i = i+3;
	}

	return longestLen;
}

int ESTAssembly::readEstQualFile() {
	vector<string> ests; //store all the ests.
	vector<string> scores; //store all the scores corresponding the the ests.
	int longestLen = 0;

	ifstream in;

	//read est file
	in.open(estFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);
		while (!in.eof()) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str.size()==0) getline(in, str);
			// first line is comment line which begins from '>'
			if (str[0] == '>') {	//comment line begins from '>'
				//vector<string> paras = split(str, '_');
				//ests.push_back(paras[1]);
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
							int theSize = estStr.size();
							if (theSize > longestLen)
								longestLen = theSize;
							break;
						}
					}
					getline(in, str);
				}
				if (str.size() == 0) {
					int theSize = estStr.size();
					if (theSize > longestLen)
						longestLen = theSize;
					ests.push_back(estStr);
				}
			}
		}
	}
	in.close();

	ifstream in2;
	//read quality file
	in2.open(qualFileName.c_str(), ios::in);

	if (in2.good()) {
		string str;
		getline(in2, str);

		while (!in2.eof()) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str.size()==0) getline(in2, str);
			// first line is comment line which begins from '>'
			if (str[0] == '>') {	//comment line begins from '>'
				//get est in the next lines
				getline(in2, str);
				string estStr;
				while (str.size() != 0) {
					str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
					if (str.size() != 0)	{
						if (str[0] != '>') {
							estStr.append(str+" ");
						} else  {
							scores.push_back(estStr);
							break;
						}
					}
					getline(in2, str);
				}
				if (str.size() == 0) {
					scores.push_back(estStr);
				}
			}
		}
	}
	in2.close();

	//generate a graph from the input files
	int i=0;
	while (i<ests.size()) {
		//the sequence is upper-case
		g->addNode(Node(ests[i], ests[i+1], toUpperCase(ests[i+2]), getIntScores(scores[i/3])));
		i = i+3;
	}
	return longestLen;
}

//given a record from a quality file, transform it to int scores. The scores in the input string are separated by one space.
std::vector<int> ESTAssembly::getIntScores(const std::string& str) {
	vector<int> ret;
	string tmp = "";
	for (int i=0; i<str.size(); i++) {
		if (str[i] != ' ') {
			tmp.push_back(str[i]);
		} else {
			ret.push_back(atoi(tmp.c_str()));
			tmp = "";
		}
	}
	if (tmp.size() > 0) {
		ret.push_back(atoi(tmp.c_str()));
	}

	return ret;
}

/*
 * read a minimum spanning tree from the input MST file
 */
void ESTAssembly::readMST() {
	int nOfNodes = g->graphNodes.size();
	if (nOfNodes == 0) {
		cout << "zero nodes in the peace-generated input MST file!" << endl;
		exit(1);
	}

	//vector<vector<int> > nodes(nOfNodes-1, vector<int>(3));	//store edges in MST, there are n-1 edges, n is number of nodes.
	vector<vector<float> > nodes;

	//read mst from the input file
	ifstream in;
	in.open(mstFileName.c_str(), ios::in);

	if (in.good()) {
		string str;
		getline(in, str);

		while (str.size() != 0) {
			str.erase(str.find_last_not_of(" \n\r\t")+1); //trim str
			if (str[0] != '#') {	//comment line begins from '#'
				vector<string> paras = split(str, ',');
				if (atoi(paras[0].c_str()) != -1) {	//-1 means root of MST
					vector<float> tmp(3);
					tmp[0] = atof(paras[0].c_str());
					tmp[1] = atof(paras[1].c_str());
					tmp[2] = atof(paras[2].c_str());
					nodes.push_back(tmp);
				}
			}
			getline(in, str);
		}
	}
	in.close();

	// Make a undirected MST.
	DefGraph mst(nOfNodes);
	for (int j=0; j<nodes.size(); j++) {
		addEdge(mst, int(nodes[j][0]), int(nodes[j][1]), nodes[j][2], false);
	}

	g->setMst(mst);
}


vector<string> ESTAssembly::split(const string& str, char delimit) {
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

string ESTAssembly::toUpperCase(const string& str) {
	string retStr;
	for (int i=0; i<str.size(); i++) {
		retStr.push_back(toupper(str[i]));
	}
	return retStr;
}
