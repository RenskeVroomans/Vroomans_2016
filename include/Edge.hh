#ifndef EdgeHeader
#define EdgeHeader

#include "Vertex.hh"

class Edge
{
 public:
  Vertex *V;
  int weight;
  double HH;
  int use;

  Edge();
  ~Edge();
  Edge(Vertex *ve,int w,int u, double H);
  Edge(const Edge &e);
};

#endif
