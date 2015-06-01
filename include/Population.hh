#ifndef PopulationHeader
#define PopulationHeader

#include "Agent.hh"
#include "Cell.hh"
#include "Network.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Header.hh"
#include "Misc.hh"
#include <cstdio>
#include <stdlib.h> 

class Population
{
 public:
  Agent* Ags[NI][NJ];
  Agent* Ags2[NI][NJ];
  int pnrags_;
  //int maxnrlonganddifbands;
  double childbins_bands[MaxNrBands+1];
  double childbins_longbands[MaxNrBands+1];
  double childbins_difbands[MaxNrBands+1];
  double childbins_nonexpfitness[2*MaxNrBands+1];
  double parentbins_bands[MaxNrBands+1];
  double parentbins_longbands[MaxNrBands+1];
  double parentbins_difbands[MaxNrBands+1];
  double parentbins_nonexpfitness[2*MaxNrBands+1];
  std::list<int> popancestry;
  typedef std::list<int>::iterator iter;


  int pagidc_;
  Population();
  ~Population();
  void InitPopulation();
  void InitPopulationFromGenome(char *fname);
  void DeathOfAgent(int i,int j);
  void ReproduceAgent(int t,int sourcei,int sourcej,int desti,int destj);
  void Margolus();
  void Shuffle();
  void Step(int t);
  void WriteFieldToFile(int time);
  void WritePopFitnessToFile(int time);
  void WriteGenomeLengthToFile(int time);
  void WriteBandsToFile(int time);
  void WriteFittestMatricesToFile(int time);
  void WriteAncestriesToFile(int t);
  void ReadAncestorsFromFile(char *fname);
  void WriteFullAncestry(Agent *C);
  void WriteAgentToFile(char *dirname,int i, int j);
};
#endif
