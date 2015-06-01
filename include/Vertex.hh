#ifndef VertexHeader
#define VertexHeader

#include "Gene.hh"

class Vertex
{
 public:
  Gene *Gen;
  int label;

  Vertex();
  ~Vertex();
  Vertex(Gene *G,int l);
  Vertex(const Vertex &v);
};

#endif
