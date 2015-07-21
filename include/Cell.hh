#ifndef CellHeader
#define CellHeader

#include "Header.hh"
#include "Genome.hh"
#include "Network.hh"
#include "Misc.hh"
#include "stdlib.h"
#include <iomanip>

class Cell
{
 public:
  int nrg;
  int age;
  int id;
  int divisioncounter;
  double genestates[300];
  double proteinstates[NrGeneTypes];//genes of same type together determine 1 protein level
  double maintproteinstates[NrGeneTypes];
  double varmaintproteinstates[NrGeneTypes];
  double signalstates[NrSignGeneTypes];

  Cell(int cid);
  ~Cell();
  // Cell& operator=(const Cell& celltocopy);
  void SetCellState(Genome *G);
  void UpdateCellState(Network *N);
  void UpdateCellAge();
};

#endif
