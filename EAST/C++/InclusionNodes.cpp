#include "InclusionNodes.h"
#include <iostream>
#include <stack>
using namespace std;

void InclusionNodes::addNode1(int chd, int parent) {
	this->nodes.insert(pair<int, int> (chd, parent));
}

void InclusionNodes::addNode2(int chd, int parent) {
	this->nodes2.insert(pair<int, int> (chd, parent));
}

int InclusionNodes::getSize() {
	return this->nodes.size();
}

bool InclusionNodes::containInclusionNode(int idx) {
	map<int, int>::iterator ret;
	ret = nodes.find(idx);
	if (ret == nodes.end()) {
		return false;
	} else {
		return true;
	}
}

vector<int> InclusionNodes::containPNode(int pIdx, int size) {
	vector<int> chdIdx;

	if (nodesGraph.size() == 0) {
		makeGraphFromAllNodes(size);
	}
	DefGraph nodesTree = Prim(nodesGraph, pIdx, true);

	stack<int> tStack;
	tStack.push(pIdx);
	while (!tStack.empty()) {
		int tmpPIdx = tStack.top();
		tStack.pop();
		for (EdgeIterator it=nodesTree[tmpPIdx].begin(); it!=nodesTree[tmpPIdx].end(); it++) {
			chdIdx.push_back((*it).node);
			tStack.push((*it).node);
		}
	}

	return chdIdx;
}

vector<int> InclusionNodes::getAllChdNodes() {
	vector<int> chdIdx;
	map<int, int>::iterator it;
	multimap<int, int>::iterator it2;
	map<int, int> indices; //Used to avoid the duplicate children index.

	for (it = nodes.begin(); it != nodes.end(); it++) {
		int idx = (*it).first;
		indices.insert(pair<int, int> (idx, 0));
	}
	for (it2 = nodes2.begin(); it2 != nodes2.end(); it2++) {
		int idx = (*it2).first;
		indices.insert(pair<int, int> (idx, 0));
	}
	for (it = indices.begin(); it != indices.end(); it++) {
		chdIdx.push_back((*it).first);
	}
	return chdIdx;
}

void InclusionNodes::printAllNodes() {
	map<int, int>::iterator it;
	multimap<int, int>::iterator it2;

	cout << "childIndex\tparentIndex";
	for (it = nodes.begin(); it != nodes.end(); it++) {
		cout << (*it).first << "\t" << (*it).second << endl;
	}
	for (it2 = nodes2.begin(); it2 != nodes2.end(); it2++) {
		cout << (*it2).first << "\t" << (*it2).second << endl;
	}
}

//Construct a graph from nodes and nodes2, store the tree into nodesGraph.
//It is called by containPNode which will return all the descendants of the input node.
void InclusionNodes::makeGraphFromAllNodes(int size) {
	for (int i=0; i<size; i++) {
		addNode(nodesGraph);
	}

	for (map<int, int>::iterator it = nodes.begin(); it != nodes.end(); it++) {
		addEdge(nodesGraph, (*it).second, (*it).first, 0, true); //parent->child, directed.
	}
	for (multimap<int, int>::iterator it2 = nodes2.begin(); it2 != nodes2.end(); it2++) {
		addEdge(nodesGraph, (*it2).second, (*it2).first, 0, true); //parent->child, directed.
	}
}
