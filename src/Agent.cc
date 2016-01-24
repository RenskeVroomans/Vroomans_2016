#include "Agent.hh"
#include "stdio.h"
#include <stdlib.h> 
#include <png.h>
#include <zlib.h>


Agent::Agent()
{
  wfspeed=1.;
}

Agent::~Agent()
{
 
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
  //if(cells!=NULL)
  //{
  cells.clear();
      // cells=NULL;
      //}
}

void Agent::CreateAgentFromScratch(int agid,int t)
{
  InitAgent(agid,0,t);
  G->GenerateGenome();
  N->BuildNetwork(G);
  //wfspeed=1.+(uniform()*0.3-0.15);//if you want to add developmental noise
  DevelopAgent();
  DetermineFitness(0);
#if defined(RUN) || defined(RERUN)
  ancestry.push_back(-1);//beginning of ancestry, no parent
#endif
}

void Agent::CreateAgentFromGenome(int agid,int t,char *fname)
{
  InitAgent(agid,0,t);
  G->CreateGenomeFromFile(fname,0);
  N->BuildNetwork(G);
  //wfspeed=1.+(uniform()*0.3-0.15);//if you want to add developmental noise
  DevelopAgent();
  DetermineFitness(0);
#if defined(RUN) || defined(RERUN)
  ancestry.push_back(-1);//beginning of ancestry, no parent
#endif
}

void Agent::CloneAgentFromAgent(Agent *A,int agid,int t)
{
  InitAgent(agid,A->agentid,t);
  G->CloneGenome(A->G);
  N->BuildNetwork(G);
  //wfspeed=1.+(uniform()*0.3-0.15);//if you want to add developmental noise
  DevelopAgent();
  DetermineFitness(0);
#if defined(RUN) || defined(RERUN)
  ancestry.push_back(-1);//beginning of ancestry, no parent
#endif
}

void Agent::CreateAgentFromParent(Agent *A,int agid,int t)
{
  InitAgent(agid,A->agentid,t);
  G->CloneGenome(A->G);
  G->MutateGenome(t);
  N->BuildNetwork(G);
  //wfspeed=1.+(uniform()*0.3-0.15);//if you want to add developmental noise
  DevelopAgent();
  DetermineFitness(0);
  //WriteGenome("testje");
#if defined(RUN) || defined(RERUN)
  ancestry=A->ancestry;//copy ancestry of parent as to be your ancestry
  ancestry.push_back(A->agentid);//add your parent to your ancestry
#endif
}


void Agent::InitAgent(int aid,int pid,int t)
{
  int i,j,k;
 
  anrcells_=0;
  agentid=aid; //agent and parent id is to trace ancestry
  parentid=pid;
  tbirth=t;
  fitness=0;

  G=NULL;
  G=new Genome();
  N=NULL;
  N=new Network();

  for(i=0;i<NrStorages;i++)//to store gene expression pattern
  {
    for(j=0;j<NrFinalCells;j++)//to determine segmentation pattern
    {
      for(k=0;k<NrGeneTypes;k++)//and its stability in last period
      {
	E[i][j][k]=0;
	if(!(j%10))
	  cellpattern[i][j/10][k]=0;
      }
    }
  }
  
  
  for(i=0;i<NrStorages;i++)//to store celltypes (derived from gene expression pattern)
  {
    for(j=0;j<NrFinalCells;j++)
    {
      types[i][j]=0;
      ages[i][j]=-1;
    }
  }
  
  //for array pattern checking!
   for(j=0;j<NrFinalCells;j++)
   {
     averseg[j]=0;
     varseg[j]=0;
   }
    
}

void Agent::UpdateCellAges()
{ 
  list<Cell>::iterator iter;
   
  for(iter=cells.begin(); iter!=cells.end();++iter)
  {
    (*iter).UpdateCellAge();
  }
  
}

void Agent::PrintCellAges()
{ 
  list<Cell>::iterator iter;
  FILE *datafile;
  int i;
  datafile=fopen("testdirs/testdata/CellAgesIds.dat", "w");
  for(iter=cells.begin(), i=0; iter!=cells.end();++iter,i++)
  {
    fprintf(datafile,"%d\t%d\n",i,(*iter).age);
  }
  fprintf(datafile,"\n");
  fclose(datafile); 
  
}
void Agent::PrintAgentData(int tijd)
{ 
  list<Cell>::iterator iter;
  char filename[60];
  snprintf(filename,60,"testdirs/testdata/agentgenes%04d.dat",tijd);
  FILE *datafile;
  datafile=fopen(filename, "w");
  int i;
  
  for(iter=cells.begin(),i=0; iter!=cells.end();++iter, i++)
  {
    fprintf(datafile, "%d\t",i);
    for(int j=0; j<NrGeneTypes;j++)
    {
      fprintf(datafile, "%.3lf\t",(*iter).proteinstates[j]);
    }
    fprintf(datafile, "\n");
  }
  
  fclose(datafile);
}

void Agent::DevelopAgent()
{
  int i,j;
  list<Cell>::iterator iter;
  FormZygote();
  
  for(i=0;i<=NrDevSteps;i++)
  {
    //printf("step %d\n",i);
    //Step 1 cell cell signalling
    CellCellSignalling(i);
    
    //Step 2 intracellular dynamics
    IntracellularDynamics();
//     #ifdef TWOMORPHS
//     if(i==NrDevSteps/2)
//     {
//       RemoveGradient();
//     }
//     #endif
    //step3 divide
    DivideCells();
    
    //Step 3 storage
    if(i%StorageInt==0){
      StoreAgentState((int)(i/StorageInt));//stores cell types and cell ages
      if(storecellspec)
	StoreCellStates((int)(i/StorageInt));
    }
    //Step 4 maintenance
    MaintenanceIntracellularDynamics(i);
  
    UpdateCellAges();
 
  }
 
}

void Agent::RemoveGradient()
{
   #ifdef TWOMORPHS
     list<Cell>::iterator iter=cells.begin();
     for (int i=0;iter!=cells.end(); ++iter,i++)
     {
       (*iter).proteinstates[1]=0;
     }
     #endif
}

void Agent::FormZygote()
{
  list<Cell>::iterator iter;
  int startnrcells=InitNrCells;
  int i;
  //setup of identical cells
  for(i=0;i<startnrcells;i++){
    Cell c(i);
    c.SetCellState(G);
    cells.push_back(c);
    anrcells_++;
  }
  
  #ifdef INITPOINT
  //Set maternal gene in posterior-most cell to 100
  iter=cells.begin();
  for (i=0;i<NrMatGeneTypes; i++)
  {
    (*iter).proteinstates[i]=100.;
  }
  #endif
  
  #ifdef INITGRAD
  iter=cells.begin();
  for (i=0;iter!=cells.end(); ++iter,i++)
  {
    (*iter).proteinstates[0]=100.*exp(-morphdecay*i);
  }
  #endif
  
  #ifdef TWOMORPHS
  iter=cells.begin();
  (*iter).proteinstates[0]=100.;
  list<Cell>::reverse_iterator riter=cells.rbegin();
  for (i=0;riter!=cells.rend(); ++riter,i++)
  {
    (*riter).proteinstates[1]=100.*exp(-0.2*i);
  }
  #endif
   
   for(iter=cells.begin();iter!=cells.end();++iter)
   {
     for(int ii=0;ii<G->gnrgenes_;ii++)
     {
       if(G->genetypeorder[ii]<NrMatGeneTypes)
	 (*iter).genestates[ii]=(*iter).proteinstates[G->genetypeorder[ii]]/(double)G->genetypenrs[G->genetypeorder[ii]];
     }  
   }
}

void Agent::CellCellSignalling(int t)
{
  int i,j, k;
  list<Cell>::iterator iter,forw,back;
  //superimposed signalling
#ifdef GRADIENT
  //production,decay and diffusion of gradient. May still require the multiple intermediate steps (using Forward Euler).
  /*  iter=cells.begin();
      for(j=0;j<NrMatGeneTypes;j++){
      for(i=0;i<growzonesize;iter++,i++)
     // (*iter).proteinstates[j]+=HT*(morphprod-Decay*(*iter).proteinstates[j]);
      MorphogenDiffusion(j);
    }*/

  //an intracellular gradient of maternal product slowly being degraded
  iter=cells.begin();
  for(i=0;iter!=cells.end();++iter,i++)
  {
    for(j=0;j<1;j++) //j<NrMatGeneTypes. Gene nr 1 for TWOMORPHS does not decay in a normal way
      {
	if(i>0){//make growzonesize be one.
	 (*iter).proteinstates[j]+=RungeKutta4(0., morphdecay,(*iter).proteinstates[j]);//-=HT*morphdecay*(*iter).proteinstates[j];
	}
	else
	  (*iter).proteinstates[j]=100.;
      }
      
      for(int ii=0;ii<G->gnrgenes_;ii++)
      {
	if(G->genetypeorder[ii]<NrMatGeneTypes)
	  (*iter).genestates[ii]=(*iter).proteinstates[G->genetypeorder[ii]]/(double)G->genetypenrs[G->genetypeorder[ii]];
      }    
  }
  
#endif

#ifdef FREEMORPH
//the morphogen behaves like every other gene; just forms asymmetric initial condition.
// iter=cells.begin();
// 
// (*iter).proteinstates[0]+=morphinflux;
// for(int ii=0;ii<G->gnrgenes_;ii++)
// {
//   if(G->genetypeorder[ii]<NrMatGeneTypes)
//     (*iter).genestates[ii]=(*iter).proteinstates[G->genetypeorder[ii]]/(double)G->genetypenrs[G->genetypeorder[ii]];
// }
#endif

#ifdef FREEDIFF
//the morphogen is regulated like other genes, but can also diffuse. to be combined with freemorph to create influx in posterior 
//and regulation by network (where also decay happens)
double temp=0;
double temparray[NrFinalCells];

for(k=0; k<Nrdiffsteps; k++)
{
  ///empty storage
  for(i=0; i<NrFinalCells; i++)
    temparray[i]=0.;
  
  //go through cells, find new values
  for(i=0,iter=cells.begin();iter!=cells.end();++iter,i++) 
  {
    temp=0;
    ///diffusion:
    forw=iter;
    ++forw;
    back=iter;
    --back;
    
    if(iter!=cells.begin())
      temp=(*back).proteinstates[0]-(*iter).proteinstates[0];//+neighbour, -self
    if(forw!=cells.end())
      temp+=(*forw).proteinstates[0]-(*iter).proteinstates[0];
    
    temparray[i]=(HT/(double)Nrdiffsteps)*(DifCoef*temp);//-morphdecay*(*iter).proteinstates[j]);
    //printf("i=%d, temparray=%.2lf\n",i,temparray[i]);
    
  }
  ///update the proteinarray
  for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
  {
    (*iter).proteinstates[0]+=temparray[i];
  }
}
///update gene states? if diffusion, perhaps unnecessary.. NO: proteinstates set to 0 in Network!!!
for(int ii=0;ii<G->gnrgenes_;ii++)
{
  if(G->genetypeorder[ii]<NrMatGeneTypes)
    (*iter).genestates[ii]=(*iter).proteinstates[G->genetypeorder[ii]]/(double)G->genetypenrs[G->genetypeorder[ii]];
}
#endif

#ifdef MORPHUP
//hard-wired upregulation of morphogen by itself
double morphstate, Hstate;
Hstate=H*H;
double mactiv=40;
for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
{
  morphstate=(*iter).proteinstates[0]*(*iter).proteinstates[0];
  (*iter).proteinstates[0]+=mactiv*morphstate/(Hstate+morphstate);
  for(int ii=0;ii<G->gnrgenes_;ii++)
  {
    if(G->genetypeorder[ii]<NrMatGeneTypes)
      (*iter).genestates[ii]=(*iter).proteinstates[G->genetypeorder[ii]]/(double)G->genetypenrs[G->genetypeorder[ii]];
  }
  
}


#endif

#ifdef POSTERIORSIGNAL
//Posterior-most cell stays on, rest off. Gene is not regulated by network
iter=cells.begin();
for(i=0;iter!=cells.end();++iter,i++) 
  for(j=0;j<NrMatGeneTypes;j++)
  {
    if(i>0)
    {
      (*iter).proteinstates[j]=0.;
    }
    else
      (*iter).proteinstates[j]=Emax;
  }
#endif

#ifdef MORPHDIFF
//the morphogen is high in the posterior, diffuses through the tissue and decays
double temp=0;
double temparray[NrFinalCells];

for(k=0; k<Nrdiffsteps; k++)
{
  //empty storage
  for(i=0; i<NrFinalCells; i++)
    temparray[i]=0.;
  
  //go through cells, find new values
  for(i=0,iter=cells.begin();iter!=cells.end();++iter,i++) 
  {
    temp=0;
    ///diffusion:
    forw=iter;
    ++forw;
    back=iter;
    --back;
    
    if(iter!=cells.begin())
      temp=(*back).proteinstates[0]-(*iter).proteinstates[0];//+neighbour, -self
    if(forw!=cells.end())
      temp+=(*forw).proteinstates[0]-(*iter).proteinstates[0];
    
    temparray[i]=(HT/(double)Nrdiffsteps)*(DifCoef*temp);//-morphdecay*(*iter).proteinstates[j]);
    //printf("i=%d, temparray=%.2lf\n",i,temparray[i]);
    
  }
  //update the array
  for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
  {
    if(i>0)
      (*iter).proteinstates[0]+=temparray[i];
    else{
      (*iter).proteinstates[0]=Emax;
      //printf("high conc\n");
    }
  }
}

///decay
iter=cells.begin();
for(i=0;iter!=cells.end();++iter,i++)
  for(j=0;j<NrMatGeneTypes;j++)
  {
    if(i>0){//make growzonesize be one.
	 (*iter).proteinstates[j]-=HT*morphdecay*(*iter).proteinstates[j];
    }
    else
      (*iter).proteinstates[j]=Emax;
  }


#endif




#ifdef WAVEFRONT  
  for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
    { gauss
      if(i>(int)((wfspeed*HT*t)))
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    (*iter).proteinstates[j]=Emax*(((double)(i-wfspeed*HT*t)*(double)(i-wfspeed*HT*t))/(((double)(i-wfspeed*HT*t)*(double)(i-wfspeed*HT*t))+(10.0*10.0)));
	}
      else
	{
	  for(j=0;j<NrMatGeneTypes;j++)
	    (*iter).proteinstates[j]=0;
	}
    }
#endif


  //evolvable signalling
  if(anrcells_>1 && NrSignGeneTypes)//only relevant for multiple cells and if any gene types are capable of doing this
    {      
      //direct cell cell signalling: notch-delta like
      //take sum of inputs from direct neighbours, no input from yourself
      //initialisation
      for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
	for(j=0;j<NrSignGeneTypes;j++)
	  (*iter).signalstates[j]=0;
      //update
      for(iter=cells.begin(),i=0;iter!=cells.end();++iter,i++)
	{
	  forw=iter;
	  back=iter;
	  forw++;
	  back--;
	  for(j=0;j<NrSignGeneTypes;j++)
	    {
	      if(i>0 && forw!=cells.end())
		(*iter).signalstates[j]=(*back).proteinstates[j+NrMatGeneTypes]+(*forw).proteinstates[j+NrMatGeneTypes];
	      else if(i==0)
		(*iter).signalstates[j]=(*forw).proteinstates[j+NrMatGeneTypes];
	      else if(forw==cells.end())//iter is at the last element
		(*iter).signalstates[j]=(*back).proteinstates[j+NrMatGeneTypes];
	    } 
	}
    }
}


void Agent::IntracellularDynamics()
{
  list<Cell>::iterator iter;
  iter=cells.begin();

  while(iter!=cells.end())
    {
      (*iter).UpdateCellState(N);
      iter++;
    }
}


void Agent::DivideCells(void)//the first cell in the growthzone divides, inserting a cell after it
{
  list<Cell>::iterator iter;
    
  
  for(iter=cells.begin(); iter!=cells.end(); ++iter)
  {
  
    #ifdef STATHEAD
    if((*iter).id==0 || (*iter).id>=InitNrCells) //the head zone does not divide
    {
      #endif
      ///animal too big: do not divide.
      if (cells.size()>=NrFinalCells) 
	return;
    
      ///reset division counter if the growth gene is below the activation threshold
	if((*iter).proteinstates[GrowGeneNr]<ThOff)
	{
	  (*iter).divisioncounter=0;
	}
	
	///if the cell has a high concentration of growth gene, increase division counter or divide.
	else if((*iter).proteinstates[GrowGeneNr]>ThOn) //>=Emax/(Decay*2)
	{
	  if(uniform()<divprob)  //(*iter).divisioncounter>=divinterval
	  {
	    (*iter).divisioncounter=0; //reset division counter
	    Cell c(anrcells_);
	    anrcells_++;
	    
	    //halve the concentration of the grow gene
	    (*iter).proteinstates[GrowGeneNr]*=0.5;
	    for(int k=0;k<G->gnrgenes_;k++)
	    {
	      if(G->genetypeorder[k]==GrowGeneNr)
	      {
		(*iter).genestates[k]*=0.5;
	      }
	    }
	    
	    //set daughter cell to state of parent cell
	    for(int i=0;i<NrGeneTypes;i++) 
	    {
	      c.proteinstates[i]=(*iter).proteinstates[i];
	      c.maintproteinstates[i]=(*iter).maintproteinstates[i];
	      c.varmaintproteinstates[i]=(*iter).varmaintproteinstates[i];
	    }
	    for(int k=0;k<G->gnrgenes_;k++)
	    {
	      c.genestates[k]=(*iter).genestates[k];
	    }
	    
	    cells.insert(iter,c);//cell is inserted before the current cell.
	    //nrofdivs++;
	  }
	  
	  else
	    (*iter).divisioncounter++;
	
	}
	#ifdef STATHEAD
    }
#endif
  }
}

void Agent::StoreCellStates(int I) //this function stores the temporal info per cell, ipv per position as StoreAgentState does.
{
  list<Cell>::iterator iter;
  iter=cells.begin();

  if(I>=NrStorages)
    {
      printf("storage array not large enough\n");
      exit(1);
    }

  while(iter!=cells.end())
    {
      if(!((*iter).id%10))
      {
	for(int k=0;k<NrGeneTypes;k++)
	cellpattern[I][(*iter).id/10][k]=(int)(*iter).proteinstates[k]; 
	
      }
      iter++;
    
    }
}
/*
void Agent::StoreAgentState(int I)
{
  //store cell states
  int j,k;
  list<Cell>::iterator iter;
  iter=cells.begin();

  if(I>=NrStorages)
    {
      printf("storage array not large enough\n");
      exit(1);
    }

  j=0;
  while(iter!=cells.end())
    {
      for(k=0;k<NrGeneTypes;k++)
	E[I][j][k]=(int)(*iter).proteinstates[k]; 

      types[I][j]=ProteinStateToCellType((*iter).proteinstates);
      ages[I][j]=(*iter).age;

      iter++;
      j++;
    }
}*/

//this one stores with fixed anterior position, ipv fixed posterior position.
void Agent::StoreAgentState(int I)
{
  //store cell states
  int j,k;
  list<Cell>::reverse_iterator riter;
  riter=cells.rbegin();

  if(I>=NrStorages)
    {
      printf("storage array not large enough\n");
      exit(1);
    }

  j=0;
  while(riter!=cells.rend())
    {
      for(k=0;k<NrGeneTypes;k++)
	E[I][j][k]=(int)(*riter).proteinstates[k]; 

      types[I][j]=ProteinStateToCellType((*riter).proteinstates);
      ages[I][j]=(*riter).age;

      ++riter;
      j++;
    }
}

void Agent::MaintenanceIntracellularDynamics(int i)
{
  int k;
  list<Cell>::iterator iter;
  //store gene expression at time=endtimedevelopment-maintenanceinterval 
  if(!averagepattern && i==NrDevSteps-MaintInt)
  { 
    for(iter=cells.begin();iter!=cells.end();++iter) //for now, consider NrFinalCells the max nr of cells.
	for(k=0;k<NrGeneTypes;k++)
	  (*iter).maintproteinstates[k]=(*iter).proteinstates[k];//maintproteinstates now contains the average.
  }
  else if(averagepattern && i>NrDevSteps-MaintInt)//-1
  {
    for(iter=cells.begin();iter!=cells.end();++iter) //for now, consider NrFinalCells the max nr of cells.
	for(k=0;k<NrGeneTypes;k++){
	  (*iter).maintproteinstates[k]+=(*iter).proteinstates[k]/(double)MaintInt;//maintproteinstates now contains the average.
	  (*iter).varmaintproteinstates[k]+=(*iter).proteinstates[k]*(*iter).proteinstates[k];
	}
  }
  if(i==NrDevSteps-MaintInt) //store how big the animal was at the start of the stability check window
    maintsize=cells.size();
  
}

/*
//this one performs the check on the position rather than within the cell
void Agent::MaintenanceIntracellularDynamics(int i)
{
  int k;
  list<Cell>::reverse_iterator riter;
  //store gene expression at time=endtimedevelopment-maintenanceinterval 
  if(!averagepattern && i==NrDevSteps-MaintInt-1)
  { 
    for(k=0, riter=cells.rbegin();riter!=cells.rend();++riter,k++) //for now, consider NrFinalCells the max nr of cells.
	averseg[k]=(*riter).proteinstates[SegmGeneNr];
  }
  else if(averagepattern && i>=NrDevSteps-MaintInt)
  {
    for(k=0, riter=cells.rbegin();riter!=cells.rend();++riter, k++) //for now, consider NrFinalCells the max nr of cells.
    {
      averseg[k]+=(*riter).proteinstates[SegmGeneNr]/(double)MaintInt;//maintproteinstates now contains the average.
      varseg[k]+=(*riter).proteinstates[SegmGeneNr]*(*riter).proteinstates[SegmGeneNr];
    }
  }
  
  if(i==NrDevSteps-MaintInt) //store how big the animal was at the start of the stability check window
    maintsize=cells.size();
  
}*/


void Agent::DetermineFitness(int mode)
{
  int i,j,k,ii, jj;
  int nrboundaries=0;
  nrbands=0;
  nrlongbands=0;
  int boundaries[NrFinalCells+1];//Stores positions of boundaries?
  list<Cell>::iterator iter, forw,ext;
  
  if((G->CheckViability())!=TRUE && !mode)
  {
    fitness=-1;
  }
  else
  {  
    /// //count the boundaries of segments (where segmentation gene goes up or down)
    //initialize the counting arrays
    nrboundaries=0;
    nrbands=0;
    nrlongbands=0;
    for(i=0;i<NrFinalCells+1;i++)//!!
       boundaries[i]=0;
    for(i=0;i<NrFinalCells;i++){
      for(j=0;j<4;j++)
	bands[i][j]=0;
    }
    k=SegmGeneNr; //k is the gene type which represents the segments (usually 5)
    
    boundaries[nrboundaries]=0;//begin of animal
    nrboundaries++; 
    
    //initialize the iterator: do we include the growth zone in the fitness yes or no?    
    iter=cells.begin();
    jj=0;
    
    //choose whether we check over an average or compare two time points
    if(averagepattern){
      for(i=jj;iter!=cells.end()&& i<(cells.size()-4);++iter,i++){
	forw=iter;
	ext=iter;
	for(j=i+1;j<i+5;j++)//boundary need not be sharp so do not jCreateGenomeust check nearest neighbour grid point
	 {
	   forw++;
	   if(forw==cells.end())
	     break;
	   else
	   {
	     if(((*iter).maintproteinstates[k]<ThOff && (*forw).maintproteinstates[k]>=ThOn) || 
	       ((*iter).maintproteinstates[k]>=ThOn && (*forw).maintproteinstates[k]<ThOff)) //there's a boundary between these two positions
	     {
	       if((j-i)==1)//boundary is right next to current position
		     boundaries[nrboundaries]=j;
	       else //try to get closer to real position
		     {
		       ii=i+1;
		       ext++;
		       while(ii<j)
		       {
			 if(((*ext).maintproteinstates[k]<ThOff && (*forw).maintproteinstates[k]>=ThOn) || 
			   ((*ext).maintproteinstates[k]>=ThOn && (*forw).maintproteinstates[k]<ThOff))  
			 {
			   iter=ext;
			   i=ii;
			 }
			 else
			   ii=1000;//probably to break out of while loop :S
			 ii++;
			 ext++;
		       }
		       if((j-i)==1)
			 boundaries[nrboundaries]=j;
		       else if((j-i)==2 || (j-i)==3)
			 boundaries[nrboundaries]=j-1;
		       else if((j-i)==4)
			 boundaries[nrboundaries]=j-2;
		     }
		   i=j-1;//boundaries[nrboundaries];
		   j=1000;//to break out of J for loop (I think). 
		   nrboundaries++;  	      
		 }
	     }
	 }
      }//cell loop
      }//if statement
      
      else 
      {
	for(i=jj;iter!=cells.end()&& i<(cells.size()-4);++iter,i++){
	  forw=iter;
	  ext=iter;
	  for(j=i+1;j<i+5;j++)//boundary need not be sharp so do not just check nearest neighbour grid point
	  {
	    forw++;
	    if(forw==cells.end())
	      break;
	    else
	    {
	      if(((*iter).proteinstates[k]<ThOff && (*forw).proteinstates[k]>=ThOn) || 
		((*iter).proteinstates[k]>=ThOn && (*forw).proteinstates[k]<ThOff)) //there's a boundary between these two positions
	      {
		if((j-i)==1)//boundary is right next to current position
		     boundaries[nrboundaries]=j;
		else //try to get closer to real position 
		{
		  ii=i+1;
		  ext++;
		  while(ii<j)
		  {
		    if(((*ext).proteinstates[k]<ThOff && (*forw).proteinstates[k]>=ThOn) || 
		      ((*ext).proteinstates[k]>=ThOn && (*forw).proteinstates[k]<ThOff))  
		     {
		      iter=ext;
		      i=ii;
		    }
		    else
		      ii=1000;//probably to break out of while loop :S
		      ii++;
		    ext++;
		  }
		  if((j-i)==1)
		    boundaries[nrboundaries]=j;
		  else if((j-i)==2 || (j-i)==3)
		    boundaries[nrboundaries]=j-1;
		  else if((j-i)==4)
		    boundaries[nrboundaries]=j-2;
		}
		i=j-1;//boundaries[nrboundaries];
		j=1000;//to break out of J for loop (I think). 
		nrboundaries++;  	      
	      }
	    }
	  }
	}//cell loop
      }

     boundaries[nrboundaries]=cells.size();//end of animal //NrFinalCells
     nrboundaries++;
     
     //Determine from the boundaries the number of bands
     //Determine from the boundaries the length of bands
     //Determine which bands are long enough to count
     int length,minlength;
     vector <int> lengtharray;
     nrbands=nrboundaries-1;
     length=0;
     minlength=minbandsize;
     for(i=0;i<nrbands;i++)
     {
       bands[i][0]=boundaries[i];//start segment
       bands[i][1]=boundaries[i+1];//end segment
       length=boundaries[i+1]-boundaries[i];
       bands[i][2]=length;
       if(length>=minlength)
       {
	 bands[i][3]=1;
// 	 /// for regularity 
// 	 if(boundaries[i]<anrcells_-InitNrCells)
// 	   lengtharray.push_back(length); //don't count the headbands
// 	 
	   
	 nrlongbands++;
       }
     }
     
     ///regularity check
     /*int median;
     int difflengths=0;
     if(lengtharray.size()>2) //need at least 3 segments to give penalty for regularity
     {
       sort(lengtharray.begin(),lengtharray.end());
       median=lengtharray[lengtharray.size()/2];
       for(i=0; i<lengtharray.size(); i++)
       {
	 difflengths+=(lengtharray[i]-median)*(lengtharray[i]-median);
       }
     }*/
     
//      double difflengths=0.;
//      if (nrlongbands>2)//do not punish the 2-stripe state for not being equal
//      {
//        for(i=0; i<nrbands-2; i++) //note that the head segment may give some issues
//        {
// 	 difflengths+=(double)abs(lengtharray[i]-lengtharray[i+1]);// /((lengtharray[i]+lengtharray[i+1])/2);
//        }
//        //difflengths/=nrlongbands;
//      }
//      else 
//        difflengths=0;
//      
     //min nr of (long enough) segments is 1
     if(nrlongbands==0)
       nrlongbands=1;
     
    
     for(i=0;i<nrboundaries-1;i++)
       bandtypes[boundaries[i]]=10;
     
     ///check stability of the pattern
     int differentcells=0;
     k=SegmGeneNr;
     iter=cells.begin();
    
     if(!averagepattern)
     {
       //compute difference between stored profile at time=enddevelopment-maintenanceinterval
       //versus e`xpression profile at time=enddevelopment
       
       for(;iter!=cells.end();++iter)
       {
	 if(((*iter).proteinstates[k]<ThOff && (*iter).maintproteinstates[k]>ThOn)||
	   ((*iter).proteinstates[k]>ThOn && (*iter).maintproteinstates[k]<ThOff))
	   differentcells++;
       }
       instpenalty=0.1*differentcells;
     }
     else{
       for(;iter!=cells.end();++iter)
       {
	 (*iter).varmaintproteinstates[k]=sqrt((*iter).varmaintproteinstates[k]/(double)MaintInt-(*iter).maintproteinstates[k]*(*iter).maintproteinstates[k]);//calculate the standard deviation at this position
	 if((*iter).varmaintproteinstates[k]>5.00)
	   differentcells++;
	 instpenalty=0.1*differentcells;
       }
     }
     
     ///various penalties
     glpenalty=genepen*G->gnrgenes_+tfbspen*G->gnrtfbs_; //0.00002 and 0.000002 for sims6_gsizepen, 0.00001 and 0.000001 normal,0.0001 and 0.00001*G for sims6_gsizepen10, 0.001 and 0.0001*G for sims6_gsizepen100  
     
     shortsegpenalty=nrbands-nrlongbands;
     
     sizefit=0.;
     if(sizebonus>0.0001)
       sizefit=sizebonus*(min((double)cells.size()-InitNrCells,(double)targetsize)-InitNrCells);//growing bigger helps by itself
     if(sizepen>0.0001)
       sizefit-=sizepen*(max(0., (double)cells.size()-(double)targetsize));//penalty for growing too large
     if(stablesizepen>0.0001)
       sizefit-=stablesizepen*(cells.size()-maintsize);//penalty for growing at the end of development; 
         
     //regularity penalty is now a bonus
//      if(nrlongbands>2)
//        regpenalty=(regpen*2)/(difflengths+2); //penalty for different sizes of the "big enough" segments
//      else
//        regpenalty=0;
     
     nonexpfitness=max(0.001,(nrlongbands+sizefit-glpenalty-instpenalty-shortsegpenalty/*+regpenalty*/));
     double selcoef=1.0;
     fitness=exp(selcoef*nonexpfitness)-1.;
     
   }
   
   //PrintCellAges();
}
/*
void Agent::DetermineFitness(int mode)
{
  int i,j,k,ii, jj;
  int nrboundaries=0;
  nrbands=0;
  nrlongbands=0;
  int boundaries[NrFinalCells+1];//Stores positions of boundaries?
  list<Cell>::reverse_iterator riter, forw,ext;
  
  if((G->CheckViability())!=TRUE && !mode)
    fitness=-1;
  else
  {  
    /// //count the boundaries of segments (where segmentation gene goes up or down)
    //initialize the counting arrays
    nrboundaries=0;
    nrbands=0;
    nrlongbands=0;
    for(i=0;i<NrFinalCells+1;i++)//!!
       boundaries[i]=0;
    for(i=0;i<NrFinalCells;i++){
      for(j=0;j<4;j++)
	bands[i][j]=0;
    }
    k=SegmGeneNr; //k is the gene type which represents the segments (usually 5)
    
    boundaries[nrboundaries]=0;//begin of animal
    nrboundaries++; 
    
    //initialize the iterator: do we include the growth zone in the fitness yes or no?    
    riter=cells.rbegin();
    jj=0;
    
    //choose whether we check over an average or compare two time points
    if(averagepattern){
      for(i=0; i<(cells.size()-4); i++) //posloop
      {
	for(j=i+1;j<i+5;j++) //neighbour loop
	{
	  if(j>=cells.size())
	    break;
	  else
	  {
	   if((averseg[i]<ThOff && averseg[j]>=ThOn) || (averseg[i]>=ThOn && averseg[j]<ThOff)) //boundary
	   {
	     if((j-i)==1)//boundary is right next to current position
	       boundaries[nrboundaries]=j;
	     else //try to get closer to real position
	     {
	       ii=i+1;
	       while(ii<j)
	       {
		 if((averseg[ii]<ThOff && averseg[j]>=ThOn) || (averseg[ii]>=ThOn && averseg[j]<ThOff))
		   i=ii;
		 else
		   break;
		 ii++;
	       }
	       if((j-i)==1)
		 boundaries[nrboundaries]=j;
	       else if((j-i)==2 || (j-i)==3)
		 boundaries[nrboundaries]=j-1;
	       else if((j-i)==4)
		 boundaries[nrboundaries]=j-2;
	     }
	     i=j-1;//boundaries[nrboundaries];
	     j=1000;//to break out of J for loop (I think). 
	     nrboundaries++;  
	   }
	  }
	}
      }
    }
    
    else 
    {
      for(i=jj;riter!=cells.rend()&& i<(cells.size()-4);++riter,i++){
	forw=riter;
	ext=riter;
	for(j=i+1;j<i+5;j++)//boundary need not be sharp so do not just check nearest neighbour grid point
	  {
	    forw++;
	    if(forw==cells.rend())
	      break;
	    else
	    {
	      if(((*riter).proteinstates[k]<ThOff && (*forw).proteinstates[k]>=ThOn) || 
		((*riter).proteinstates[k]>=ThOn && (*forw).proteinstates[k]<ThOff)) //there's a boundary between these two positions
	      {
		if((j-i)==1)//boundary is right next to current position
		     boundaries[nrboundaries]=j;
		else //try to get closer to real position 
		{
		  ii=i+1;
		  ext++;
		  while(ii<j)
		  {
		    if(((*ext).proteinstates[k]<ThOff && (*forw).proteinstates[k]>=ThOn) || 
		      ((*ext).proteinstates[k]>=ThOn && (*forw).proteinstates[k]<ThOff))  
		     {
		      riter=ext;
		      i=ii;
		    }
		    else
		      ii=1000;//probably to break out of while loop :S
		      ii++;
		    ext++;
		  }
		  if((j-i)==1)
		    boundaries[nrboundaries]=j;
		  else if((j-i)==2 || (j-i)==3)
		    boundaries[nrboundaries]=j-1;
		  else if((j-i)==4)
		    boundaries[nrboundaries]=j-2;
		}
		i=j-1;//boundaries[nrboundaries];
		j=1000;//to break out of J for loop (I think). 
		nrboundaries++;  	      
	      }
	    }
	  }
	}//cell loop
      }

     boundaries[nrboundaries]=cells.size();//end of animal //NrFinalCells
     nrboundaries++;
     
     //Determine from the boundaries the number of bands
     //Determine from the boundaries the length of bands
     //Determine which bands are long enough to count
     int length,minlength;
     vector <int> lengtharray;
     nrbands=nrboundaries-1;
     length=0;
     minlength=minbandsize;
     for(i=0;i<nrbands;i++)
     {
       bands[i][0]=boundaries[i];//start segment
       bands[i][1]=boundaries[i+1];//end segment
       length=boundaries[i+1]-boundaries[i];
       bands[i][2]=length;
       if(length>=minlength)
       {
	 bands[i][3]=1;
	 /// for regularity 
	 if(boundaries[i]<anrcells_-InitNrCells)
	   lengtharray.push_back(length); //don't count the headbands
	 
	   
	 nrlongbands++;
       }
     }
     
     ///regularity check
     int median;
     int difflengths=0;
     if(lengtharray.size()>2) //need at least 3 segments to give penalty for regularity
     {
       sort(lengtharray.begin(),lengtharray.end());
       median=lengtharray[lengtharray.size()/2];
       for(i=0; i<lengtharray.size(); i++)
       {
	 difflengths+=(lengtharray[i]-median)*(lengtharray[i]-median);
       }
     }

     
     //min nr of (long enough) segments is 1
     if(nrlongbands==0)
       nrlongbands=1;
     
    
     for(i=0;i<nrboundaries-1;i++)
       bandtypes[boundaries[i]]=10;
     
     ///check stability of the pattern
     int differentcells=0;
     k=SegmGeneNr;

     riter=cells.rbegin();
    
     if(!averagepattern)
     {
       //compute difference between stored profile at time=enddevelopment-maintenanceinterval
       //versus e`xpression profile at time=enddevelopment
       
       for(i=0;riter!=cells.rend();++riter,i++)
       {
	 if(((*riter).proteinstates[k]<ThOff && averseg[i]>ThOn)||
	   ((*riter).proteinstates[k]>ThOn && averseg[i]<ThOff))
	   differentcells++;
       }
       instpenalty=0.1*differentcells;
     }
     else{
       for(i=0;i<cells.size();i++)
       {
	 varseg[i]=sqrt(varseg[i]/(double)MaintInt-averseg[i]*averseg[i]);//calculate the standard deviation at this position
	 if(varseg[i]>5.00)
	   differentcells++;
	 instpenalty=0.1*differentcells;
       }
     }
     
     ///various penalties
     glpenalty=genepen*G->gnrgenes_+tfbspen*G->gnrtfbs_; //0.00002 and 0.000002 for sims6_gsizepen, 0.00001 and 0.000001 normal,0.0001 and 0.00001*G for sims6_gsizepen10, 0.001 and 0.0001*G for sims6_gsizepen100  
     
     shortsegpenalty=nrbands-nrlongbands;
     
     sizefit=sizebonus*(min((double)cells.size()-InitNrCells,(double)targetsize)-InitNrCells)-sizepen*(max(0., (double)cells.size()-(double)targetsize))-stablesizepen*(cells.size()-maintsize); //growing bigger helps by itself
          
     regpenalty=regpen*difflengths; //penalty for different sizes of the "big enough" segments
     
     nonexpfitness=max(0.001,(nrlongbands+sizefit-glpenalty-instpenalty-shortsegpenalty-regpenalty));
     double selcoef=1.0;
     fitness=exp(selcoef*nonexpfitness)-1.;
     
   }
   
   //PrintCellAges();
}*/



void Agent::WriteGenome(char *dirname)
{
  FILE *f;
  char fname[800];
  Genome::iter it;
  Gene *gene;
  TFBS *tfbs;
  
  sprintf(fname,"%s/%s/CodedGenomeAgent%.10d",despath,dirname,agentid);
  //sprintf(fname,"%s/CodedGenomeAgent%.10d",dirname,agentid);
  f=fopen(fname,"w");
  for(it=G->ChromBBList->begin();it!=G->ChromBBList->end();it++)
    {
      if(G->IsGene(*it))
	{
	  gene=dynamic_cast<Gene *>(*it);
	  fprintf(f,"%i ",0);//to indicate that it is a gene
	  fprintf(f,"%i ",(gene->type));//type of gene
	  fprintf(f,"%f ",(gene->DD));//decay rate of gene
	  fprintf(f,"%f ",(gene->EE));//decay rate of gene
	}
      else if(G->IsTFBS(*it))
	{
	  tfbs=dynamic_cast<TFBS *>(*it);
	  if(tfbs->weight<0)
	    fprintf(f,"%i ",-1);//to indicate that it is a repressing tfbs
	  else
	    fprintf(f,"%i ",1);//to indicate that it is a activating tfbs
	  fprintf(f,"%i ",(tfbs->type));//type of tfbs
	  fprintf(f,"%f ",(tfbs->HH));//Hill constant of tfbs
	}
    }
  fclose(f);
}



int Agent::CellTypeToColor(int type)
{
  int color;
  //  int modulo=251;//prime

  if(type==0)
    color=254;
  else
    color=(type*40)%251;
  return color;
}

int Agent::SegmentToColor(int type,int max)
{
  int color;
  //  int modulo=251;//prime

  color=(int)((255./(double)max)*(double)type);
  
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
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*(NrStorages);//NrStorages+1
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  // int type;
  int color;

  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      {
	color=CellTypeToColor(types[i][j]);
	for(ii=0;ii<zoom;ii++)
	  for(jj=0;jj<zoom;jj++)
	    celltypes[zoom*i+ii][zoom*j+jj]=color;
      }

//   //extra row to show segment generalized celltype classes
//   i=NrStorages;
//   for(j=0;j<NrFinalCells;j++)
//     {
//       color=CellTypeToColor(bandtypes[j]);
//       for(ii=0;ii<zoom;ii++)
// 	for(jj=0;jj<zoom;jj++)
// 	  celltypes[zoom*(NrStorages)+ii][zoom*j+jj]=color;
//     }
 
  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
      {
	RGBdata[i*WW*3+j*3+0]=colormap[celltypes[i][j]][0];
	RGBdata[i*WW*3+j*3+1]=colormap[celltypes[i][j]][1];
	RGBdata[i*WW*3+j*3+2]=colormap[celltypes[i][j]][2];
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
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//   // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);
}

void Agent::WriteDivisionProfile(char *dirname, int c)
{
  int i,j;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=4;
  const int WW=zoom*NrFinalCells;
  const int LL=zoom*(NrStorages);
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  // int type;
  int color;

  for(i=0;i<NrStorages;i++)
    for(j=0;j<NrFinalCells;j++)
      {
	color=ages[i][j];
	for(ii=0;ii<zoom;ii++)
	  for(jj=0;jj<zoom;jj++)
	    celltypes[zoom*i+ii][zoom*j+jj]=color;
      }

   
  for(i=0;i<LL;i++)
    for(j=0;j<WW;j++)
    {
     /* if(celltypes[i][j]>=0) //graded coloring of cell age
      {
	RGBdata[i*WW*3+j*3+0]=max(0.0, 255.-celltypes[i][j]);
	RGBdata[i*WW*3+j*3+1]=max(0.0, 255.-2*celltypes[i][j]);
	RGBdata[i*WW*3+j*3+2]=max(0.0, 255.-4*celltypes[i][j]);
      }*/
      if(celltypes[i][j]>-1 && celltypes[i][j]<StorageInt) 
      {
	RGBdata[i*WW*3+j*3+0]=255;
	RGBdata[i*WW*3+j*3+1]=255;
	RGBdata[i*WW*3+j*3+2]=255;
      }
      else if(celltypes[i][j]>=StorageInt)
      {
	RGBdata[i*WW*3+j*3+0]=0;
	RGBdata[i*WW*3+j*3+1]=0;
	RGBdata[i*WW*3+j*3+2]=0;
      }
      else //non-tissue part
      {
	RGBdata[i*WW*3+j*3+0]=127;
	RGBdata[i*WW*3+j*3+1]=127;
	RGBdata[i*WW*3+j*3+2]=127;
	
      }
    }
    
  if(!c)  
    sprintf(fname,"%s/%s/GrowthAgent%.10d.png",despath,dirname,agentid);
  else
    sprintf(fname,"%s/%s/GrowthAgent%.10d_%d.png",despath,dirname,agentid,c);
  
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
  png_set_IHDR(png_ptr, info_ptr,WW,LL,
	       8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
	       PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
//   // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++ ) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
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

 /*
  k=0;
  for(k=0; k<NrStorages; k++)
  {
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
  }*/
  
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
 
  ////Temporal Dynamics///
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
  
  for(j=0;j<NrFinalCells/10;j++)
  {
    sprintf(fname,"%s/%s/TemporalDynamics%.10d_CellId%i",despath,dirname,agentid,j*10);
    f=fopen(fname,"w");
    for(i=0;i<NrStorages;i++)
    {
      fprintf(f,"%i\t",i);
      for(k=0;k<NrGeneTypes;k++)
	fprintf(f,"%i\t",cellpattern[i][j][k]);
      fprintf(f,"\n");
    }
    fclose(f);
  }
  
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
  fprintf(f,"bodysize %i\n",cells.size());
  fprintf(f,"nonexpfitness %f\n",nonexpfitness);
  fprintf(f,"fitness %f\n",fitness);
  fprintf(f,"glpenalty %f\n",glpenalty);
  fprintf(f,"instpenalty %f\n",instpenalty);
  fprintf(f,"shortsegpenalty %f\n",shortsegpenalty);
  fprintf(f,"sizefitness %f\n",sizefit);
  
  for(i=0;i<nrbands;i++)
    {
      j=(int)(bands[i][0]+0.5*bands[i][2]);
      fprintf(f,"band %i begin %i end %i length %i long enough %i\t type %i\n",i,bands[i][0],bands[i][1],bands[i][2],bands[i][3],bandtypes[j]);
    }
  fclose(f);
}
