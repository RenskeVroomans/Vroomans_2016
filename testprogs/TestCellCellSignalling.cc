#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <list>

using namespace std;

class Cell
{

public:
  Cell(int ident):id(ident){}
  int id;
  double signalstates;
  double proteinstates;
};


int main()
{
  list<Cell> cellijst;
  list<Cell>::iterator iter, forw, back;
  list<Cell>::reverse_iterator riter;
  int i=0;
  
  
  
  /**initialize cellijst**/
  while(i<10){
    cellijst.push_back(Cell(i));
    riter=cellijst.rbegin();
    (*riter).proteinstates=(double)(i)*0.1;
    i++;
  }
  
 // for(iter=cellijst.begin();iter!=cellijst.end();++iter)
  //  cout << "cell "<<(*iter).id << ", concentration: "<<(*iter).proteinstates<<"\n"<<endl;
  
    
  /** do integration **/
    
  int NrDifSteps=100;
  int NrSteps=1000;
  double HT=0.01;
  double DifCoef=0.5;
  
  for(int l=0; l<NrSteps; l++){
    for(int k=0;k<NrDifSteps;k++)
    {
      //initialisation
      for(iter=cellijst.begin(),i=0;iter!=cellijst.end();++iter,i++)
	(*iter).signalstates=0;
      //update
	for(iter=cellijst.begin(),i=0;iter!=cellijst.end();++iter,i++)
	{
	  forw=iter;
	  back=iter;
	  forw++;
	  back--;
	  
	  if(i>0)//input from left cell
	  {
	    (*iter).signalstates+=(*back).proteinstates;//+neighbor
	    (*iter).signalstates-=(*iter).proteinstates;//-self
	  }
	  if(forw!=cellijst.end())//input from right cell 
	  {
	    (*iter).signalstates+=(*forw).proteinstates;//+neighbor
	    (*iter).signalstates-=(*iter).proteinstates;//-self
	  }
	  (*iter).signalstates=HT*DifCoef*(*iter).signalstates;//this is only the change in concentration
	}
	//compute new conc, store in proteinstates
	for(iter=cellijst.begin(),i=0;iter!=cellijst.end();++iter,i++)
	  (*iter).proteinstates=(*iter).proteinstates+(*iter).signalstates; //change is added here
    }  
    //assign new conc to signalling
    for(iter=cellijst.begin(),i=0;iter!=cellijst.end();++iter,i++)
      (*iter).signalstates=(*iter).proteinstates;//the signalling state is the new protein concentration in that cell
      
      
      for(iter=cellijst.begin();iter!=cellijst.end();++iter)
	cout <<(*iter).id << "\t"<<(*iter).proteinstates<<endl;
      cout <<"\n\n";
  }
  
  return 0;
}