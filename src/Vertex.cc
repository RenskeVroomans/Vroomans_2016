#include "Vertex.hh"
#include "Header.hh"

Vertex::Vertex()
{
  Gen=NULL;
  label=-1;
}


Vertex::~Vertex()
{
}

Vertex::Vertex(Gene *G,int l)
{
  Gen=G;
  label=l;
}

Vertex::Vertex(const Vertex &v)
{
  Gen=v.Gen;
  label=v.label;
}
