#include "Prim.h"
#include <math.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <assert.h>
#include <limits.h>

using namespace std;

void addNode(DefGraph& G) {
  G.push_back(vector<Edge>());
}

void addEdge(DefGraph& G, int s, int t, float weight, int reverse, bool directed) {
  assert(s < (int)G.size() && t < (int)G.size());
  G[s].push_back(Edge(t,weight,reverse));
  if (!directed)
    G[t].push_back(Edge(s,weight,reverse));
}

int degree(DefGraph G, int n) {
  assert(n < (int)G.size());
  return G[n].size();
}

