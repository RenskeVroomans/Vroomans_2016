#ifndef EdgeHeader
#define EdgeHeader

#include "Vertex.hh"

class Edge
{
 public:
  Vertex *V;
  int weight;
  int use;

  Edge();
  ~Edge();
  Edge(Vertex *ve,int w,int u);
  Edge(const Edge &e);
};

#endif
