#include "Cell.hh"
#include "stdio.h"


Cell::Cell(int cid)
{
  int i;
  age=0;
  id=cid;
  divisioncounter=0;
  for(i=0;i<NrGeneTypes;i++)
  {
    proteinstates[i]=0;
    maintproteinstates[i]=0;
    varmaintproteinstates[i]=0;
  }
  for(i=0;i<NrSignGeneTypes;i++)
    signalstates[i]=0;
  
}

Cell::~Cell()
{
}

void Cell::SetCellState()
{
  int i;
  
 
  for(i=0;i<NrGeneTypes;i++)
    {
      if(i<NrMatGeneTypes)
	proteinstates[i]=0;
      //else if(i>SegmGeneNr &&i<=SegmGeneNr+3)//switch on 3 genes after the segm gene
      else if (i==6 || i==7)//switch on the signalling genes and the growth gene //NrMatGeneTypes+NrSignGeneTypes+1
	proteinstates[i]=100.0;
      else
	proteinstates[i]=0.0;
    }
}
 
void Cell::UpdateCellState(Network *N)
{
  N->UpdateNetworkState(1,proteinstates,signalstates);
}

void Cell::UpdateCellAge()
{
  age++;
}



