#ifndef NetworkHeader
#define NetworkHeader

#include "Header.hh"
#include "Genome.hh"
#include "ChromBB.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Vertex.hh"
#include "Edge.hh"
#include "Misc.hh"


class Network
{
 public:
  list<Vertex *> *VL;
  typedef std::list<Vertex *>::iterator iterv;
  list<list <Edge *> > *AL;
  typedef list<list <Edge *> >::iterator iterel;
  typedef list<Edge *>::iterator itere;
  typedef std::list<int>::iterator iternsf;
  typedef std::list<int>::reverse_iterator iternsb;
 

  Network();
  ~Network();
  void BuildNetwork(Genome *G);
  void UpdateNetworkState(int steps,double genestates[300],double proteinstates[NrGeneTypes],double signalstates[NrSignGeneTypes]);
 };

#endif
