#include "Edge.hh"
#include "Header.hh"

Edge::Edge()
{
  V=NULL;
  weight=0;
  use=0;
}


Edge::~Edge()
{
}

Edge::Edge(Vertex *ve,int w, int u)
{
  V=ve;
  weight=w;
  use=u;
}

Edge::Edge(const Edge &e)
{
  V=e.V;
  weight=e.weight;
}
