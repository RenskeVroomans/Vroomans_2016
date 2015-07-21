#include "Edge.hh"
#include "Header.hh"

Edge::Edge()
{
  V=NULL;
  weight=0;
  use=0;
  HH=0.;
}


Edge::~Edge()
{
}

Edge::Edge(Vertex *ve,int w, int u, double H)
{
  V=ve;
  weight=w;
  use=u;
  HH=H;
}

Edge::Edge(const Edge &e)
{
  V=e.V;
  HH=e.HH;
  weight=e.weight;
}
