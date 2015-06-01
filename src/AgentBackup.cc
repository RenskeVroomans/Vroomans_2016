#include "Agent.hh"
#include "stdio.h"
#include <stdlib.h> 
#include <png.h>
#include <zlib.h>


Agent::Agent()
{
}

Agent::~Agent()
{
  int i;

  if(G!=NULL)
    {
      delete (G);
      G=NULL;
    }
  if(N!=NULL)
    {
      delete (N);
      N=NULL;
    }
  for(i=0;i<NrFinalCells;i++)
    {
      if(CL[i]!=NULL)
	{
	  delete CL[i];
	  CL[i]=NULL;	
	}
    }  
}

void Agent::CreateAgentFromScratch(int agid,int t)
{
  InitAgent(agid,0,t);
  G->GenerateGenome();
  N->BuildNetwork(G);
  DevelopAgent();
  DetermineFitness();
#ifdef RUN
  ancestry.push_back(-1);//beginning of ancestry, no parent
#endif
}

void Agent::CloneAgentFromAgent(Agent *A,int agid,int t)
{
  InitAgent(agid,A->agentid,t);
  G->CloneGenome(A->G);
  N->BuildNetwork(G);
  DevelopAgent();
  DetermineFitness();
#ifdef RUN
  ancestry.push_back(-1);//beginning of ancestry, no parent
#endif
}

void Agent::CreateAgentFromParent(Agent *A,int agid,int t)
{
  InitAgent(agid,A->agentid,t);
  G->CloneGenome(A->G);
  G->MutateGenome(t);
  N->BuildNetwork(G);
  DevelopAgent();
  DetermineFitness();
#ifdef RUN
  ancestry=A->ancestry;//copy ancestry of parent as to be your ancestry
  ancestry.push_back(A->agentid);//add your parent to your ancestry
#endif
}


void Agent::InitAgent(int aid,int pid,int t)
{
  int i,j,k;
 
  anrcells_=0;
  agentid=aid;
  parentid=pid;
  tbirth=t;
  fitness=0;

  G=NULL;
  G=new Genome();
  N=NULL;
  N=new Network();

  for(i=0;i<NrFinalCells;i++)
    CL[i]=NULL;
  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      for(k=0;k<NrGeneTypes;k++)
	E[i][j][k]=0;
  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      types[i][j]=0;
}

void Agent::DevelopAgent()
{
  int i;

  FormZygote();
  for(i=0;i<=NrDevSteps;i++)
    {
      //printf("time %i nrcels %i\n",i,anrcells_);
      //Step 1 cell cell signalling
      CellCellSignalling(i);
      //Step 2 intracellular dynamics
      IntracellularDynamics();
      //Step 3 countdown cell division tracker
      CountDownDivisionTracker();
      //Step 4 growth
      GrowthByDivision();
      //Step 5 storage
      if(i%StorageInt==0)
	StoreAgentState((int)(i/StorageInt));
      //Step 6 averaging
      AverageIntracellularDynamics(i);
      //Step 7 maintenance
      MaintenanceIntracellularDynamics(i);
    }
}

void Agent::FormZygote()
{
#ifdef SINGLECELL
  CL[0]=new Cell();
  CL[0]->SetCellState();
  anrcells_++;
  CL[0]->divisioncounter=MinDivTime;//(int)(uniform()*MaxDevDivTime)+MinDivTime;
#endif
#ifdef MULTIPLECELLS
  int startnrcells=InitNrCells;//10;
  int i,j;
  //setup of identical cells
  for(i=0;i<startnrcells;i++)
    {
      CL[i]=new Cell();
      CL[i]->SetCellState();
      anrcells_++;
      CL[i]->divisioncounter=MinDivTime;//(int)(uniform()* MaxDevDivTime)+MinDivTime;
    }
  //induce asymmetry
#ifdef INTERNALASYMMETRY
  for(i=0;i<startnrcells;i++)
    {
      if(i<startnrcells/8)
	CL[i]->proteinstates[NrMatGeneTypes+NrSignGeneTypes+1]=100;
      else
	CL[i]->proteinstates[NrMatGeneTypes+NrSignGeneTypes+1]=0;
    }
#endif
#ifdef MATERNALASYMMETRY
#ifdef BLOCK
   for(i=0;i<startnrcells;i++)
     {
       if(i<startnrcells/8)
	 {
	   for(j=0;j<NrMatGeneTypes;j++)
	     CL[i]->proteinstates[j]=100;
	 }
       else
	 {
	    for(j=0;j<NrMatGeneTypes;j++)
	      CL[i]->proteinstates[j]=0;
	 }
     }
#endif
#ifdef GRADIENT
   for(i=0;i<startnrcells;i++)
    for(j=0;j<NrMatGeneTypes;j++)
      CL[i]->proteinstates[j]=Emax*exp(-i/20.0);//Emax*exp(-i/16.0);
#endif
#ifdef WAVEFRONT
   for(i=0;i<startnrcells;i++)
     for(j=0;j<NrMatGeneTypes;j++)
       CL[i]->proteinstates[j]=Emax*(((double)i*(double)i)/(((double)i*(double)i)+(10.0*10.0)));
#endif
#ifdef REVERSEWAVEFRONT
   for(i=0;i<startnrcells;i++)
     for(j=0;j<NrMatGeneTypes;j++)
       CL[i]->proteinstates[j]=Emax*((10.0*10.0)/(((double)i*(double)i)+(10.0*10.0)));
#endif
#endif

#endif
 
}


void Agent::CellCellSignalling(int t)
{
  int i,j;
#ifdef MULTIPLECELLS
#ifdef MATERNALASYMMETRY
#ifdef BLOCK
  for(i=0;i<anrcells_;i++)
    {
      if(i<anrcells_/8)
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    CL[i]->proteinstates[j]=__gnu_cxx::power(0.925,t)*100;
	}
      else
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    CL[i]->proteinstates[j]=0;
	}
    }
#endif
#ifdef GRADIENT
  for(i=0;i<anrcells_;i++)
    for(j=0;j<NrMatGeneTypes;j++)
      CL[i]->proteinstates[j]=__gnu_cxx::power(0.925,t)*Emax*exp(-i/20.0);//__gnu_cxx::power(0.925,t)*Emax*exp(-i/16.0);
#endif
#ifdef WAVEFRONT  
  for(i=0;i<anrcells_;i++)
    {
      if(i>(int)((HT*t)))
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    CL[i]->proteinstates[j]=Emax*(((double)(i-HT*t)*(double)(i-HT*t))/(((double)(i-HT*t)*(double)(i-HT*t))+(10.0*10.0)));
	}
      else
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    CL[i]->proteinstates[j]=0;
	}
    }
#endif
#ifdef REVERSEWAVEFRONT
   for(i=0;i<anrcells_;i++)
     {
       if(i<(int)((HT*t)))
	 {
	   for(j=0;j<NrMatGeneTypes;j++)
	     CL[i]->proteinstates[j]=Emax;
	 }
       else
	 {
	   for(j=0;j<NrMatGeneTypes;j++)
	     CL[i]->proteinstates[j]=Emax*((10.0*10.0)/(((double)(i-HT*t)*(double)(i-HT*t))+(10.0*10.0)));
	 }
     }
#endif

#endif
#endif

#ifdef DIRECTSIGNALLING
  if(anrcells_>1)
    {
      for(i=0;i<anrcells_;i++)
	for(j=0;j<NrSignGeneTypes;j++)
	  CL[i]->signalstates[j]=0;
      for(i=0;i<anrcells_;i++)
	for(j=0;j<NrSignGeneTypes;j++)
	  {
	    if(i>0 && i<anrcells_-1)
	      CL[i]->signalstates[j]=CL[i-1]->proteinstates[j+NrMatGeneTypes]+CL[i+1]->proteinstates[j+NrMatGeneTypes];
	    else if(i==0)
	      CL[i]->signalstates[j]=CL[i+1]->proteinstates[j+NrMatGeneTypes];
	    else if(i==anrcells_-1)
	      CL[i]->signalstates[j]=CL[i-1]->proteinstates[j+NrMatGeneTypes];
	  } 
    }
#endif
#ifdef DIFFUSIVESIGNALLING
  if(anrcells_>1)
    {      
      double DifCoef=0.3;//0.1;      
      int NrDifSteps=1;//10;//50:
      for(j=0;j<NrSignGeneTypes;j++)
	{
	  for(int k=0;k<NrDifSteps;k++)
	    {
	      //initialisation
	      for(i=0;i<anrcells_;i++)
		CL[i]->signalstates[j]=0;
	      //update
	      for(i=0;i<anrcells_;i++)
		{
		  if(i>0)//input from left cell
		    {
		      CL[i]->signalstates[j]+=CL[i-1]->proteinstates[j+NrMatGeneTypes];//+neighbor
		      CL[i]->signalstates[j]-=CL[i]->proteinstates[j+NrMatGeneTypes];//-self
		    }
		  if(i<anrcells_-1)//input from right cell
		    {
		      CL[i]->signalstates[j]+=CL[i+1]->proteinstates[j+NrMatGeneTypes];//+neighbor
		      CL[i]->signalstates[j]-=CL[i]->proteinstates[j+NrMatGeneTypes];//-self
		    }
		  CL[i]->signalstates[j]=HT*DifCoef*CL[i]->signalstates[j];
		}
	      //compute new conc, store in proteinstates
	      for(i=0;i<anrcells_;i++)
		CL[i]->proteinstates[j+NrMatGeneTypes]=CL[i]->proteinstates[j+NrMatGeneTypes]+CL[i]->signalstates[j];
	    }  
	  //assign new conc to signalling
	  for(i=0;i<anrcells_;i++)
	    CL[i]->signalstates[j]=CL[i]->proteinstates[j+NrMatGeneTypes];
	}
    }
#endif
#ifdef MIXEDSIGNALLING
   if(anrcells_>1)
    {      
      //direct cell cell signalling: notch-delta
      //take sum of inputs from direct neighbours
      //initialisation
      for(i=0;i<anrcells_;i++)
	for(j=0;j<NrSignGeneTypes/2;j++)
	  CL[i]->signalstates[j]=0;
      //update
      for(i=0;i<anrcells_;i++)
	for(j=0;j<NrSignGeneTypes/2;j++)
	  {
	    if(i>0 && i<anrcells_-1)
	      CL[i]->signalstates[j]=CL[i-1]->proteinstates[j+NrMatGeneTypes]+CL[i+1]->proteinstates[j+NrMatGeneTypes];
	    else if(i==0)
	      CL[i]->signalstates[j]=CL[i+1]->proteinstates[j+NrMatGeneTypes];
	    else if(i==anrcells_-1)
	      CL[i]->signalstates[j]=CL[i-1]->proteinstates[j+NrMatGeneTypes];
	  } 
     
      //cell signalling gene 1 /second half:
      //diffusion based cell cell signalling: Wnt, TFG-beta
      double DifCoef=0.3;//0.1;      
      int NrDifSteps=1;//10;//50:

      for(j=NrSignGeneTypes/2;j<NrSignGeneTypes;j++)
	{
	  for(int k=0;k<NrDifSteps;k++)
	    {
	      //initialisation
	      for(i=0;i<anrcells_;i++)
		CL[i]->signalstates[j]=0;
	      //update
	      for(i=0;i<anrcells_;i++)
		{
		  if(i>0)//input from left cell
		    {
		      CL[i]->signalstates[j]+=CL[i-1]->proteinstates[j+NrMatGeneTypes];//+neighbor
		      CL[i]->signalstates[j]-=CL[i]->proteinstates[j+NrMatGeneTypes];//-self
		    }
		  if(i<anrcells_-1)//input from right cell
		    {
		      CL[i]->signalstates[j]+=CL[i+1]->proteinstates[j+NrMatGeneTypes];//+neighbor
		      CL[i]->signalstates[j]-=CL[i]->proteinstates[j+NrMatGeneTypes];//-self
		    }
		  CL[i]->signalstates[j]=HT*DifCoef*CL[i]->signalstates[j];
		}
	      //compute new conc, store in proteinstates
	      for(i=0;i<anrcells_;i++)
		CL[i]->proteinstates[j+NrMatGeneTypes]=CL[i]->proteinstates[j+NrMatGeneTypes]+CL[i]->signalstates[j];
	    }  
	  //assign new conc to signalling
	  for(i=0;i<anrcells_;i++)
	    CL[i]->signalstates[j]=CL[i]->proteinstates[j+NrMatGeneTypes];
	}
    }
#endif
 
}


void Agent::IntracellularDynamics()
{
  int j;
  for(j=0;j<anrcells_;j++)
    {
      CL[j]->UpdateCellState(N);//updates genestates, proteinstates, but not signalstates
    }
}


void Agent::CountDownDivisionTracker()
{
  int j;
  for(j=0;j<anrcells_;j++)
    {
      if(CL[j]->divisioncounter>0)
	CL[j]->divisioncounter--;
    }
}

void Agent::GrowthByDivision()
{
  int j;
  int divnr;
  for(j=0;j<anrcells_;j++)
    {
      if(CL[j]->divisioncounter==0 && anrcells_<NrFinalCells)
	{
	  if(anrcells_==1)
	    divnr=1;//first asymmetric division
	  else
	    divnr=2;//other symmetric divisions
	  CellDivision(j,divnr);
	  CL[j]->divisioncounter=MinDivTime;//(int)(uniform()* MaxDevDivTime)+MinDivTime;
	  j++;//new cell not allowed to divide immediately
	  CL[j]->divisioncounter=MinDivTime;//(int)(uniform()* MaxDevDivTime)+MinDivTime;
	}
    }
}


void Agent::CellDivision(int cellnr,int divnr)
{
  int i;
  double a;
  double nrmol;
  double mu;
  double sigma;
  double dx;
  int ix;
  int ixx;

//#define DIVISIONNOISE

  //move cells to the right of dividing cell with one position to create space for new cell
  for(i=anrcells_-1;i>=cellnr+1;i--)//do this from back to front to not loose info!
    CL[i+1]=CL[i];
  //add 1 cell to nr of cells in organism
  anrcells_++;
  //create new cell to right of to be divided cell
  CL[cellnr+1]=new Cell();
  //divide protein concentrations over dividing cells
  a=1.0;//cleavage division//a=0.5;//normal division
  for(i=0;i<NrGeneTypes;i++)
    {
      if(divnr==1 && i==(NrMatGeneTypes+NrSignGeneTypes+1))
	{             
	  CL[cellnr+1]->proteinstates[i]=CL[cellnr]->proteinstates[i];
	  CL[cellnr]->proteinstates[i]=0;
	}
      else 
	{
#ifdef DIVISIONNOISE
	  nrmol=CL[cellnr]->proteinstates[i];
	  mu=nrmol*a;
	  sigma=sqrt(nrmol)*0.5;
	  dx=gauss[(int)(uniform()*10000)]*sigma+mu;
	  ix=(int)(dx+0.5);
	  ixx=mu+(mu-ix);
	  if(ix<0)
	    ix=0;
	  if(ixx<0)
	    ixx=0;
	  CL[cellnr+1]->proteinstates[i]=ix;
	  CL[cellnr]->proteinstates[i]=ixx;
#else
	  CL[cellnr+1]->proteinstates[i]=a*CL[cellnr]->proteinstates[i];//(int)(a*CL[cellnr]->proteinstates[i]+0.5);
	  CL[cellnr]->proteinstates[i]=a*CL[cellnr]->proteinstates[i];//(int)(a*CL[cellnr]->proteinstates[i]+0.5);
#endif
	}
    }
}



void Agent::StoreAgentState(int I)
{
  //store cell states
  int j,k;
  
  if(I>=NrStorages)
    {
      printf("storage array not large enough\n");
      exit(1);
    }

  for(j=0;j<anrcells_;j++)
    for(k=0;k<NrGeneTypes;k++)
      E[I][j][k]=(int)CL[j]->proteinstates[k]; 
  
  for(j=0;j<anrcells_;j++)
    types[I][j]=ProteinStateToCellType(CL[j]->proteinstates);

}

void Agent::AverageIntracellularDynamics(int i)
{
  int j,k;


  if(i==NrDevSteps-50)
    {
      //printf("nr cells %i\n",anrcells_);

      for(j=0;j<NrFinalCells;j++)
	for(k=0;k<NrGeneTypes;k++)
	  CL[j]->avgproteinstates[k]=0;
    }
  if(i>=NrDevSteps-50)
    {
      for(j=0;j<NrFinalCells;j++)
	for(k=0;k<NrGeneTypes;k++)
	  CL[j]->avgproteinstates[k]+=CL[j]->proteinstates[k];
    }
  if(i==NrDevSteps)
    {
      for(j=0;j<NrFinalCells;j++)
	for(k=0;k<NrGeneTypes;k++)
	  CL[j]->avgproteinstates[k]/=51;    
    }
}  


void Agent::MaintenanceIntracellularDynamics(int i)
{
  int j,k;
  if(i==NrDevSteps-MaintInt)
    {
      for(j=0;j<NrFinalCells;j++)
	for(k=0;k<NrGeneTypes;k++)
	  CL[j]->maintproteinstates[k]=CL[j]->proteinstates[k];
    }
}


int Agent::DifferentExpressionPattern(double p1[NrGeneTypes],double p2[NrGeneTypes])
{
  int different;
  int i;

  //0: the same
  //1: different
  different=0;
  for(i=NrMatGeneTypes;i<NrGeneTypes;i++)
    {
      if(G->idimpact[i]==1)//this gene type/protein contributes to cell id
	{
	  if( (p1[i]<ThOff && p2[i]>=ThOn) || (p1[i]>=ThOn && p2[i]<ThOff) )
	    {
	      different=1;
	      break;
	    }
	}
    }
  return different;
}


void Agent::DetermineFitness()
{
  //#define SEGMENTIDENTITY
  //#define SEGMENTS
#define FIRSTSEGMENTSTHANID

#ifdef SEGMENTIDENTITY
  int i,j,k,l,ii;
  int nrboundaries=0;
  nrbands=0;
  nrlongbands=0;
  nrlonganddifbands=0;
  int boundaries[NrFinalCells+1];//!!

  if((G->CheckViability())!=TRUE)
    fitness=-1;
  else
    {  
      for(i=0;i<NrFinalCells+1;i++)//!!
	boundaries[i]=0;
      for(i=0;i<NrFinalCells;i++)
	bandtypes[i]=0;
      for(i=0;i<NrFinalCells;i++)
	for(j=0;j<4;j++)
	  bands[i][j]=0;
      for(i=0;i<NrFinalCells;i++)
	for(j=0;j<NrGeneTypes;j++)
	  avgexprbands[i][j]=0;
      
      //Determine the boundaries between sufficiently different and near enough cells
      boundaries[nrboundaries]=0;//begin of animal
      nrboundaries++; 
      for(i=0;i<NrFinalCells-4;i++)
	{
	  for(j=i+1;j<i+5;j++)
	    {
	      //if(DifferentExpressionPattern(CL[i]->avgproteinstates,CL[j]->avgproteinstates)==1)
	      if(DifferentExpressionPattern(CL[i]->proteinstates,CL[j]->proteinstates)==1)
		{
		  //determine exact location of boundary
		  if((j-i)==1)
		    boundaries[nrboundaries]=j;
		  else
		    {
		      ii=i+1;
		      while(ii<j)
			{
			  //if(DifferentExpressionPattern(CL[ii]->avgproteinstates,CL[j]->avgproteinstates)==1)//!!
			  if(DifferentExpressionPattern(CL[ii]->proteinstates,CL[j]->proteinstates)==1)
			    i=ii;
			  else
			    ii=1000;
			  ii++;
			}
		      if((j-i)==1)
			boundaries[nrboundaries]=j;
		      else if((j-i)==2 || (j-i)==3)
			boundaries[nrboundaries]=j-1;
		      else if((j-i)==4)
			boundaries[nrboundaries]=j-2;
		    }
		  //no need to look further for boundary relative to point i,
		  //so jump out of j loop, go to new i which should be the
		  //new boundary point, the -1 is to compensate for the i++
		  i=boundaries[nrboundaries]-1;//!!
		  j=1000;
		  nrboundaries++;  	      
		}
	    }
	}
      boundaries[nrboundaries]=NrFinalCells;//end of animal 
      nrboundaries++;
 
      //Determine from the boundaries the number of bands
      //Determine from the boundaries the length of bands
      //Determine which bands are long enough to count
      //Determine the avg expression pattern in a band
      nrbands=nrboundaries-1;
      int length=0;
      int minlength=minbandsize;
      for(i=0;i<nrbands;i++)
	{
	  bands[i][0]=boundaries[i];//start segment
	  bands[i][1]=boundaries[i+1];//end segment
	  length=boundaries[i+1]-boundaries[i];
	  bands[i][2]=length;
	  if(length>=minlength)
	    {
	      bands[i][3]=1;
              nrlongbands++;
	    }
	  for(k=0;k<NrGeneTypes;k++)
	    {
	      for(j=boundaries[i];j<boundaries[i+1];j++)
		avgexprbands[i][k]+=CL[j]->proteinstates[k];//CL[j]->avgproteinstates[k];//!!
	      avgexprbands[i][k]/=length;
	    }
	}
      //determine of the long enough bands how many have 
      //a distinct expression pattern
      int first=0;
      int samefound;
      double dummy;
      int samebands[NrFinalCells];
      for(i=0;i<NrFinalCells;i++)
	samebands[i]=-1;

      for(i=0;i<nrbands;i++)
	{
	  if(bands[i][3]==1)//long enough
	    {
	      if(first==0)//first long enough band
		{
		  samebands[i]=i;
		  nrlonganddifbands++;
		  first=1;		
	        }
	      else//other long enough bands
	 	{
		  samefound=0;
		  for(j=i-1;j>=0;j--)//compare to all previous bands
		    {
		      if(bands[j][3]==1)//those that are also long enough
			{
			  if(DifferentExpressionPattern(avgexprbands[i],avgexprbands[j])==1)//?problem??
			    ;//difference found
			  else//no difference found
			    {
			      samefound=1;//bands i and j are thus the same
			      samebands[i]=samebands[j];//band i is thus of same type as previous band j
			      //because you check all previous bands you will end up with the same number
			      //as the earliest band of that type
			    }
			}
		    }			
		  if(samefound==0)//no band of same type found
		      {
			samebands[i]=i;//in that case band i first of its kind
			nrlonganddifbands++;
		      }
		}
	    } 
	}

      for(i=0;i<nrbands;i++)
	{
	  for(j=bands[i][0];j<bands[i][1];j++)
	    bandtypes[j]=samebands[i];
	}
      for(i=0;i<nrboundaries-1;i++)
	 {
	   bandtypes[boundaries[i]]=10;
	 }


      int difference=0;
      int differentcells=0;
      for(i=0;i<NrFinalCells;i++)
	{
	  for(k=0;k<NrGeneTypes;k++)
	    difference+=abs((int)CL[i]->proteinstates[k]-(int)CL[i]->maintproteinstates[k]);
	  if(difference>10)
	    differentcells++;
	    
	}
      instpenalty=0.01*differentcells;

      glpenalty=0.00001*G->gnrgenes_+0.000001*G->gnrtfbs_;

      shortsegpenalty=nrbands-nrlongbands;

      nonexpfitness=max(0.001,(sqrt(nrlonganddifbands*nrlongbands)-shortsegpenalty-glpenalty-instpenalty));
      double selcoef=1.0;//5.0;//1.0;
      fitness=exp(selcoef*nonexpfitness)-1;
    }
#endif
#ifdef SEGMENTS
  int i,j,k,ii;  if((G->CheckViability())!=TRUE)
    fitness=-1;
  else
    {
      int nrboundaries;
      int boundaries[NrFinalCells+1];
      k=SegmGeneNr;
      nrboundaries=0;
      nrbands=0;
      nrlongbands=0;
      nrlonganddifbands=0;
      for(i=0;i<NrFinalCells;i++)
	boundaries[i]=0;
      for(i=0;i<NrFinalCells;i++)
	for(j=0;j<4;j++)
	  bands[i][j]=0;
      for(i=0;i<NrFinalCells;i++)
	bandtypes[i]=0;
      
      boundaries[nrboundaries]=0;//begin of animal
      nrboundaries++; 
      for(i=0;i<NrFinalCells-4;i++)
	for(j=i+1;j<i+5;j++)
	  {
	    //if((CL[i]->avgproteinstates[k]<ThOff && CL[j]->avgproteinstates[k]>=ThOn) || 
	    // (CL[i]->avgproteinstates[k]>=ThOn && CL[j]->avgproteinstates[k]<ThOff))
	    if((CL[i]->proteinstates[k]<ThOff && CL[j]->proteinstates[k]>=ThOn) || 
	       (CL[i]->proteinstates[k]>=ThOn && CL[j]->proteinstates[k]<ThOff))  
	      {
		if((j-i)==1)
		  boundaries[nrboundaries]=j;
		else
		  {
		    ii=i+1;
		    while(ii<j)
		      {
			//if((CL[ii]->avgproteinstates[k]<ThOff && CL[j]->avgproteinstates[k]>=ThOn) || 
			// (CL[ii]->avgproteinstates[k]>=ThOn && CL[j]->avgproteinstates[k]<ThOff))
			if((CL[ii]->proteinstates[k]<ThOff && CL[j]->proteinstates[k]>=ThOn) || 
			   (CL[ii]->proteinstates[k]>=ThOn && CL[j]->proteinstates[k]<ThOff))  
			  i=ii;
			else
			  ii=1000;
			ii++;
		      }
		    if((j-i)==1)
		      boundaries[nrboundaries]=j;
		    else if((j-i)==2 || (j-i)==3)
		      boundaries[nrboundaries]=j-1;
		    else if((j-i)==4)
		      boundaries[nrboundaries]=j-2;
		  }
		i=boundaries[nrboundaries];
		j=1000;
		nrboundaries++;  	      
	      }
	  }
      boundaries[nrboundaries]=NrFinalCells;//end of animal
      nrboundaries++;
      
      //Determine from the boundaries the number of bands
      //Determine from the boundaries the length of bands
      //Determine which bands are long enough to count
      nrbands=nrboundaries-1;
      int length=0;
      int minlength=minbandsize;
      for(i=0;i<nrbands;i++)
	{
	  bands[i][0]=boundaries[i];//start segment
	  bands[i][1]=boundaries[i+1];//end segment
	  length=boundaries[i+1]-boundaries[i];
	  bands[i][2]=length;
	  if(length>=minlength)
	    {
	      bands[i][3]=1;
	      nrlongbands++;
	    }
	}

      for(i=0;i<nrboundaries-1;i++)
	 {
	   bandtypes[boundaries[i]]=10;
	 }
      
      nrlonganddifbands=nrlongbands;//we are not really looking at different bands
      
  
      int difference=0;
      int differentcells=0;
      for(i=0;i<NrFinalCells;i++)
	{
	  for(k=0;k<NrGeneTypes;k++)
	    difference+=abs((int)CL[i]->proteinstates[k]-(int)CL[i]->maintproteinstates[k]);
	  if(difference>10)
	    differentcells++;
	    
	}
      instpenalty=0.01*differentcells;
      
      glpenalty=0.00001*G->gnrgenes_+0.000001*G->gnrtfbs_;

      shortsegpenalty=nrbands-nrlongbands;


      //printf("%i %i %f %f\n",nrbands,nrlongbands,instpenalty,glpenalty);
      
      nonexpfitness=max(0.001,(nrlongbands-shortsegpenalty-glpenalty-instpenalty));
      double selcoef=1.0;//5.0;//1.0;
      fitness=exp(selcoef*nonexpfitness)-1;
    }
#endif
#ifdef FIRSTSEGMENTSTHANID
  int i,j,k,l,ii;
   
 if((G->CheckViability())!=TRUE)
    fitness=-1;
  else
    {
      int boundaries[NrFinalCells+1];
      int nrboundaries;
      k=SegmGeneNr;
      nrboundaries=0;
      nrbands=0;
      nrlongbands=0;
      nrlonganddifbands=0;
      nrhomogbands=0;
      nrdifbutrepeatedbands=0;
      for(i=0;i<NrFinalCells;i++)
	boundaries[i]=0;
      for(i=0;i<NrFinalCells;i++)
	for(j=0;j<4;j++)
	  bands[i][j]=0;
      for(i=0;i<NrFinalCells;i++)
	for(j=0;j<NrGeneTypes;j++)
	  avgexprbands[i][j]=0;
      boundaries[nrboundaries]=0;//begin of animal
      nrboundaries++; 
      for(i=0;i<NrFinalCells-4;i++)
	for(j=i+1;j<i+5;j++)
	  {
	    if((CL[i]->proteinstates[k]<ThOff && CL[j]->proteinstates[k]>=ThOn) || 
	       (CL[i]->proteinstates[k]>=ThOn && CL[j]->proteinstates[k]<ThOff))
	      {
		if((j-i)==1)
		  boundaries[nrboundaries]=j;
		else
		  {
		    ii=i+1;
		    while(ii<j)
		      {
			if((CL[ii]->proteinstates[k]<ThOff && CL[j]->proteinstates[k]>=ThOn) || 
			   (CL[ii]->proteinstates[k]>=ThOn && CL[j]->proteinstates[k]<ThOff))  
			  i=ii;
			else
			  ii=1000;
			ii++;
		      }
		    if((j-i)==1)
		      boundaries[nrboundaries]=j;
		    else if((j-i)==2 || (j-i)==3)
		      boundaries[nrboundaries]=j-1;
		    else if((j-i)==4)
		      boundaries[nrboundaries]=j-2;
		  }
		i=boundaries[nrboundaries];
		j=1000;
		nrboundaries++;  	      
	      }
	  }
      boundaries[nrboundaries]=NrFinalCells;//end of animal
      nrboundaries++;
      
      //Determine from the boundaries the number of bands
      //Determine from the boundaries the length of bands
      //Determine which bands are long enough to count
      nrbands=nrboundaries-1;
      int length=0;
      int minlength=minbandsize;
      for(i=0;i<nrbands;i++)
	{
	  bands[i][0]=boundaries[i];//start segment
	  bands[i][1]=boundaries[i+1];//end segment
	  length=boundaries[i+1]-boundaries[i];
	  bands[i][2]=length;
	  if(length>=minlength)
	    {
	      bands[i][3]=1;
	      nrlongbands++;
	    }
	  for(k=0;k<NrGeneTypes;k++)
	    {
	      for(j=boundaries[i]+1;j<boundaries[i+1]-1;j++)
		avgexprbands[i][k]+=CL[j]->proteinstates[k];//CL[j]->avgproteinstates[k];//!!
	      avgexprbands[i][k]/=(length-2);
	    }
	}
      
      //so now we know the number of bands and the number of long enough
      //bands formed by differential expression of gene k=NrSignGeneTypes+2
      //now we need to find out the number of long enough bands with a unique
      //and homogeneous expression pattern
      //1) within segment expression should be the same
      //2) between segments expression should be different
      int homog[NrFinalCells];
      for(i=0;i<NrFinalCells;i++)
	homog[i]=1;//first assume they are homogeneous than check
      for(i=0;i<nrbands;i++)//determine homogeneity of bands
	{
	  for(j=bands[i][0]+1;j<bands[i][1]-1;j++)//compare each cell in band but not on boundary
	    {
	      for(k=j+1;k<bands[i][1]-1;k++)//against all next cells in band but not on boundary
		{
		  if(DifferentExpressionPattern(CL[j]->proteinstates,CL[k]->proteinstates)==1)
		    {
		      homog[i]=0;//if difference found set to not homogeneous
		      k=1000;//jump out of k loop
		      j=1000;//jump out of j loop
		    }
		}
	    }
	}
      for(i=0;i<nrbands;i++)
	{
	  if(homog[i]==1 && bands[i][3]==1)//homogeneous and long enough
	    nrhomogbands++;
	}
      

      
      //determine of the long enough and homogeneous bands how many have 
      //a distinct expression pattern and which ones have
      //qualitatively the same expression pattern
      int first=0;
      int samefound;
      double dummy;
      int samebands[NrFinalCells];
      for(i=0;i<NrFinalCells;i++)
	samebands[i]=-1;
      
      for(i=0;i<nrbands;i++)
	{
	  if(homog[i]==1 && bands[i][3]==1)//long enough and homog
	    {
	      if(first==0)//first long enough band
		{
		  samebands[i]=i;
		  nrlonganddifbands++;
		  first=1;		
		}
	      else//other long enough bands
		{
		  samefound=0;
		  for(j=i-1;j>=0;j--)//compare to all previous bands
		    {
		      if(homog[j]==1 && bands[j][3]==1)//those that are also long enough and homog
			{
			  if(DifferentExpressionPattern(avgexprbands[i],avgexprbands[j])==1)
			    ;//difference found
			  else//no difference found
			    {
			      samefound=1;//bands i and j are thus the same
			      samebands[i]=samebands[j];//band i is thus of same type as previous band j
			      //because you check all previous band you will end up with the same number
			      //as the earliest band of that type
			    }
			}
		    }			
		  if(samefound==0)//no band of same type found
		    {
		      samebands[i]=i;//in that case band i first of its kind
		      nrlonganddifbands++;
		    }
		}
	    } 
	}
      

      //select for unique AP identity rather than different segments:
      //in this case either a segment id should occur only once, or, if it occurs
      //multiple times it should only occur in consequtive segments forming a single AP
      //region but not further along the AP axis leading to inconclusive AP identity
      //note however, that if we count all genes except maternal gene 0 and segmentation
      //gene 7 for segment identity, using this method no other genes can partake in segment
      //boundary formation, as we do not allow for differentiation within AP region but 
      //intepret this than as different AP ids that are used repeatedly
      //to solve this maybe we also should not count signalling genes and a number of random
      //genes not in AP identity
      //B

      /*nrlonganddifbands=0;
      int I;
      int nrdifbands=0;
      for(i=0;i<nrbands;i++)
	{
	  samefound=0;
	  I=i;
	  if(samebands[i]==I)//first segment with that id
	    {
	      nrdifbands++;
	      while(samebands[i]==I && i<nrbands)//find conseq segments with same id
		i++;
	      i--;
	     
	      for(j=i+1;j<nrbands;j++)//look if later non consequtive segment has same id
		{
		  if(samebands[j]==I)
		    samefound=1;
		}
	      if(samefound==0)//if not it is a unique AP region identity
		nrlonganddifbands++;
	    }
	}
      nrdifbutrepeatedbands=nrdifbands-nrlonganddifbands;
      */

      //C
      nrlonganddifbands=0;
      nrdifbutrepeatedbands=0;
      int nrdifbands=0;
      for(i=0;i<nrbands;i++)
	{
	  if(i==0)
	    {
	      if(samebands[i]==-1)
		;
	      else if(samebands[i]==i)
		nrlonganddifbands++;
	    }
	  else
	    {
	      if(samebands[i]==-1)
		;//not long enough or not homogeneous
	      else if(samebands[i]==i)
		nrlonganddifbands++;//begin of new domain
	      else if(samebands[i]==samebands[i-1])
		;//still same domain
	      else//not new, not same as before, not -1, so repetition
		{
		  nrlonganddifbands++;//rest of bands is all one domain
		  break;
		}
	    }
	}
      

      for(i=0;i<nrbands;i++)
	{
	  for(j=bands[i][0];j<bands[i][1];j++)
	    bandtypes[j]=samebands[i];
	}
      for(i=0;i<nrboundaries-1;i++)
	bandtypes[boundaries[i]]=10;


      int difference=0;
      int differentcells=0;
      for(i=0;i<NrFinalCells;i++)
	{
	  for(k=0;k<NrGeneTypes;k++)
	    difference+=abs((int)CL[i]->proteinstates[k]-(int)CL[i]->maintproteinstates[k]);
	  if(difference>4)
	    differentcells++;
	}
      instpenalty=0.1*differentcells;//0.01*differentcells;

      glpenalty=0.00001*G->gnrgenes_+0.000001*G->gnrtfbs_;

      shortsegpenalty=0;//nrbands-nrlongbands;

      //nonexpfitness=max(0.001,(nrlonganddifbands-glpenalty-instpenalty));
      nonexpfitness=max(0.001,(nrhomogbands+nrlonganddifbands-glpenalty-instpenalty));
      double selcoef=1.0;//5.0;//1.0;
      fitness=exp(selcoef*nonexpfitness)-1;

      //printf("%i %i %i %i %f %f %f %f %f %f\n",nrbands,nrlongbands,nrhomogbands,nrlonganddifbands,segmfitness,domainfitness,instpenalty,glpenalty,nonexpfitness,fitness);

    }
#endif
}



void Agent::WriteGenome(char *dirname)
{
  FILE *f;
  char fname[800];
  Genome::iter it;
  Gene *gene;
  TFBS *tfbs;
  int i;
 
  sprintf(fname,"%s/%s/CodedGenomeAgent%.10d",despath,dirname,agentid);
  f=fopen(fname,"w");
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  fprintf(f,"%i ",0);//to indicate that it is a gene
	  fprintf(f,"%i ",(gene->type));//type of gene
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight<0)
	    fprintf(f,"%i ",-1);//to indicate that it is a repressing tfbs
	  else
	    fprintf(f,"%i ",1);//to indicate that it is a activating tfbs
	  fprintf(f,"%i ",(tfbs->type));//type of tfbs
	}
    }
  for(i=0;i<NrGeneTypes;i++)
    {
      if(G->idimpact[i]==0)
	fprintf(f,"%i ",-100);
      else
	fprintf(f,"%i ",100);
    }
  fclose(f);
}



int Agent::CellTypeToColor(int type)
{
  int i;
  int color;
  int modulo=251;//prime

  if(type==0)
    color=254;
  else
    color=(type*40)%251;
  return color;
}



int Agent::ProteinStateToCellType(double proteinstates[NrGeneTypes])
{
  int i;
  int type;

  type=0;
  for(i=0;i<NrGeneTypes;i++)
    {
      if(proteinstates[i]<ThOff)
	type+=0;
      else if(proteinstates[i]<ThOn)
	type+=__gnu_cxx::power(3,i);
      else
	type+=2*__gnu_cxx::power(3,i);
    }
  return type;
}



void Agent::WriteEmbryology(char *dirname)
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*(NrStorages+1);
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  int type;
  int color;

  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      {
	color=CellTypeToColor(types[i][j]);
	for(ii=0;ii<zoom;ii++)
	  for(jj=0;jj<zoom;jj++)
	    celltypes[zoom*i+ii][zoom*j+jj]=color;
      }

  //extra row to show segment generalized celltype classes
  i=NrStorages;
  for(j=0;j<NrFinalCells;j++)
    {
      color=CellTypeToColor(bandtypes[j]);
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[zoom*(NrStorages)+ii][zoom*j+jj]=color;
    }
 
  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
      {
	RGBdata[j*LL*3+i*3+0]=colormap[celltypes[i][j]][0];
	RGBdata[j*LL*3+i*3+1]=colormap[celltypes[i][j]][1];
	RGBdata[j*LL*3+i*3+2]=colormap[celltypes[i][j]][2];
      }

  sprintf(fname,"%s/%s/CellTypesAgent%.10d.png",despath,dirname,agentid);
  PNGFileP = fopen(fname, "wb");

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,(png_voidp) NULL,
				    (png_error_ptr) NULL, 
				    (png_error_ptr) NULL );
  if(!png_ptr)
    {
      printf("out of memory\n");
      exit(1);
    }
  info_ptr = png_create_info_struct ( png_ptr );
  if(!info_ptr)
    {
      png_destroy_write_struct(&png_ptr, NULL);
      printf("out of memory\n");
      exit(1);
    }
  png_init_io ( png_ptr, PNGFileP );
  png_set_IHDR(png_ptr, info_ptr,LL,WW,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = WW-1; row >= 0; row-- ) 
    {
      row_pointer = ( RGBdata + LL * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);
}



void Agent::WriteSignalProfiles(char *dirname)
{
  
  FILE *f;
  char fname[800];
  int i,j,k;

 
  k=0;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=1;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=2;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=3;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=4;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=5;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=6;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=NrStorages/4;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=NrStorages/2;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=3*NrStorages/4;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f);  
  k=NrStorages-1;
  sprintf(fname,"%s/%s/SignalProfiles%.10d_time%i",despath,dirname,agentid,k*StorageInt);
  f=fopen(fname,"w");
  for(i=0;i<NrFinalCells;i++)
    {      
      fprintf(f,"%i\t",i);
      for(j=0;j<NrGeneTypes;j++)
	fprintf(f,"%i\t",E[k][i][j]);
      fprintf(f,"\n"); 
    }
  fclose(f); 
 

  j=0;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=10;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=20;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=30;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=40;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=50;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=60;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=70;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=80;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=90;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
  j=99;
  sprintf(fname,"%s/%s/TemporalDynamics%.10d_Cell%i",despath,dirname,agentid,j);
  f=fopen(fname,"w");
  for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",E[i][j][k]);
      fprintf(f,"\n");
    }
  fclose(f);
}


void Agent::WriteFitnessDetails(char *dirname)
{
  
  FILE *f;
  char fname[800];
  int i,j;
  
  sprintf(fname,"%s/%s/FitnessDetails%.10d",despath,dirname,agentid);
  f=fopen(fname,"w");
  fprintf(f,"nrbands %i\n",nrbands);
  fprintf(f,"nrlongbands %i\n",nrlongbands);
  fprintf(f,"nrhomogbands %i\n",nrhomogbands);
  fprintf(f,"nrlonganddifbands %i\n",nrlonganddifbands);
  fprintf(f,"nrdifbutrepeatedbands %i\n",nrdifbutrepeatedbands);
  fprintf(f,"nonexpfitness %f\n",nonexpfitness);
  fprintf(f,"fitness %f\n",fitness);
  fprintf(f,"glpenalty %f\n",glpenalty);
  fprintf(f,"instpenalty %f\n",instpenalty);
  fprintf(f,"shortsegpenalty %f\n",shortsegpenalty);
  for(i=0;i<nrbands;i++)
    {
      j=(int)(bands[i][0]+0.5*bands[i][2]);
      fprintf(f,"band %i begin %i end %i length %i long enough %i\t type %i\n",i,bands[i][0],bands[i][1],bands[i][2],bands[i][3],bandtypes[j]);
    }
  fclose(f);
}
