#include "Population.hh"


Population::Population()
{
  int i;

  for(i=0;i<MaxPopSize;i++)
    {
      Ags[i]=NULL;
      Ags2[i]=NULL;
    }
}

Population::~Population()
{
  int i;

  for(i=0;i<MaxPopSize;i++)
    {
      if((Ags[i])!=NULL)
	delete (Ags[i]);
      if((Ags2[i])!=NULL)
	delete (Ags2[i]);
    }
}


void Population::InitPopulation()
{
  int i;
  Agent *A;
  Agent *Arest;
  pagidc_=1;
  int unfit;
#define HOMOGENEOUSSTART
  //#define HETEROGENEOUSSTART
  
  for(i=0;i<MaxPopSize;i++)
    {
#ifdef HOMOGENEOUSSTART
      if(i==0)
	{
	  unfit=TRUE;
	  while(unfit==TRUE)
	    {
	      A=new Agent();
	      A->CreateAgentFromScratch(pagidc_,0);
	      if(A->fitness>0)
		{
		  Ags[i]=A;
		  pagidc_++;
		  unfit=FALSE;
		}
	      else
		delete A;
	    }
	}
      else
	{
	  Arest=new Agent();
	  Arest->CloneAgentFromAgent(A,pagidc_,0);
	  Ags[i]=Arest;
	  pagidc_++;
	}
#endif
#ifdef HETEROGENEOUSSTART
      unfit=TRUE;      
      while(unfit==TRUE)
	{
	  A=new Agent();
	  A->CreateAgentFromScratch(pagidc_,0);
	  if(A->fitness>0)
	    {
	      Ags[i]=A;
	      pagidc_++;
	      unfit=FALSE;
	    }
	  else
	    delete A;
	}
#endif
    }
}


void Population::DeathOfAgent(int i)
{
  if(Ags[i]!=NULL)
    {
      delete (Ags[i]);
      Ags[i]=NULL;
    } 
}

void Population::ReproduceAgent(int t,int sourcei,int desti)
{
  Agent *P;
  Agent *C;
 
  P=Ags[sourcei];
  if(P==NULL)
    {
      printf("parent does not exist %i\n",sourcei);
      exit(1);
    }
  //else
  //printf("location parent %i id parent %i\n",sourcei,P->agentid);
  C=new Agent();
  C->CreateAgentFromParent(P,pagidc_,t);
  Ags2[desti]=C;
  pagidc_++;
}
	  

void Population::Step(int t)
{
  int i;
  int k;
  int ii;
  double maxfit;
  double fit;
  int memi;

 
  //put Ags2 to zero so no old ones there
  for(i=0;i<MaxPopSize;i++)
    Ags2[i]=NULL;

  //make new agent for all positions in field
  for(i=0;i<MaxPopSize;i++)
    {
      //draw randomly X agents from current field and pick fittest as parent
      maxfit=-1000;
      memi=-1;
      for(k=0;k<((int)(0.25*MaxPopSize));k++)
	{
	  ii=(int)(uniform()*MaxPopSize);
	  fit=Ags[ii]->fitness;
	  if(fit>maxfit)
	    {
	      maxfit=fit;
	      memi=ii;
	    }
	}
      //let it produce child and put this child temporarily in Ags2
      ReproduceAgent(t,memi,i);
    }

  //all old agents die: remove from Ags
  for(i=0;i<MaxPopSize;i++)
    DeathOfAgent(i);

  //all new agents are replacing them: update Ags with Ags2
  for(i=0;i<MaxPopSize;i++)
    Ags[i]=Ags2[i];

}


void Population::WritePopFitnessToFile(int time)
{
  FILE *f;
  char fname[800];
  meanfit=0;
  minfit=1000;
  maxfit=-1000;
  stdfit=0;
  int i;
  int fitnessbins[150];
  int index;
 
  for(i=0;i<150;i++)
    fitnessbins[i]=0;
  
  for(i=0;i<MaxPopSize;i++)
    {
      meanfit+=Ags[i]->fitness;
      stdfit+=Ags[i]->fitness*Ags[i]->fitness;
      if(Ags[i]->fitness<minfit)
	minfit=Ags[i]->fitness;
      if(Ags[i]->fitness>maxfit)
	maxfit=Ags[i]->fitness;

      index=(int)(Ags[i]->fitness);
      if(index>=0 && index<150)
	fitnessbins[index]++;
    }
  
  meanfit/=MaxPopSize;
  stdfit/=MaxPopSize;
  stdfit=stdfit-meanfit*meanfit;
  if(stdfit>0)
    stdfit=sqrt(stdfit);
 
  sprintf(fname,"%s/PopFitnessDynamics",despath);
  f=fopen(fname,"a");
  fprintf(f,"%i\t%f\t%f\t%f\t%f\n",time,meanfit,minfit,maxfit,stdfit);
  fclose(f);

  sprintf(fname,"%s/PopFitnessHisto%.4d",despath,time);
  f=fopen(fname,"w");
  for(i=0;i<150;i++)
    fprintf(f,"%i\t%i\n",i,fitnessbins[i]);
  fclose(f);

}

void Population::WriteFitAgentToFile(int time)
{
  int i;
  for(i=0;i<MaxPopSize;i++)
    {
      if((int)(1000*Ags[i]->fitness)==(int)(1000*maxfit))
	{
	  //printf("writing agent %i with fitness %f to file\n",Ags[i]->agentid,Ags[i]->fitness);
	  //for(int j=0;j<128;j++)
	  //printf("%i ",Ags[i]->celltypes[j]);
	  //printf("\n");
	  //Ags[i]->WriteGenome();
	  Ags[i]->WriteEmbryology();
	  break;
	}
    }
}
