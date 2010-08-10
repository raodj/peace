#include <vector>

#ifndef JK_PRIM_H
#define JK_PRIM_H

class Edge {
public:
  Edge(int node, float weight, int reverse) {this->node = node; this->weight = weight; this->reverse=reverse;}
  int node;
  float weight;
  int reverse;
};


// Graph: we will use the adjacency matrix representation
//    G[i] will be the adjacency list for node i
//    G[i][j] is the jth edge touching node i.  (Order is arbitrary.)
//    G[i][j].node is the node on the other end of the edge.
//    G[i][j].distance is the weight of the edge.
typedef std::vector<std::vector<Edge> > DefGraph;
typedef std::vector<Edge>::iterator EdgeIterator;

void addNode(DefGraph& G);
void addEdge(DefGraph& G, int s, int t, float weight, int reverse, bool directed = true);
int degree(DefGraph G, int n); // "Out degree" if directed

#endif
