#ifndef GenomeHeader
#define GenomeHeader

#include "Header.hh"
#include "ChromBB.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include <typeinfo>

class Genome
{
 public:
  std::list<ChromBB*> *ChromBBList;
  typedef std::list<ChromBB*>::iterator iter;
  typedef std::list<ChromBB*>::reverse_iterator reviter;
  int glength_;
  int gnrgenes_;
  int gnrtfbs_;

  Genome();
  ~Genome();
  void CopyPartOfGenome(iter begin,iter end);
  void CopyPartOfGenomeToTempList(iter begin,iter end,list<ChromBB*> &ChromBBListTemp);
  void CloneGenome(const Genome *g);
  void GenerateGenome();
  void CreateGenomeFromFile(char *fname, int agentid); //if agentid=0, you use the function as meant in World
  bool IsGene(ChromBB *cbb) const;
  bool IsTFBS(ChromBB *cbb) const;
  int MutateGenome(int t);
  iter FindFirstTFBSInFrontOfGene(iter ii) const;
  iter FindRandomGenePosition() const;
  iter GeneMutate(iter);
  int OnlyCopy(iter ii);
  iter TFBSMutate(iter);
  void TFBSInnovation();
  int CheckViability();
  void ListContent();
  
  //for analysis
  void RemoveSegmentationGene();
  void PruneGenome();
  void DeleteGene(int nr);
  void DeleteTFBS(int nr);
};
#endif

