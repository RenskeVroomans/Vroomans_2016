#include "Population.hh"


Population::Population()
{
  int i,j;
  pnrags_=0;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	Ags[i][j]=NULL;
	Ags2[i][j]=NULL;
      }
}

Population::~Population()
{
  int i,j;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if((Ags[i][j])!=NULL)
	  delete (Ags[i][j]);
	if((Ags2[i][j])!=NULL)
	  delete (Ags2[i][j]);
      }
}


void Population::InitPopulation()
{
  int i,j;
  int ii,jj;
  Agent *A;
  Agent *Arest;
  pagidc_=1;
  int unfit;
  int counter;
 
  counter=0;
  ii=NI/2;
  jj=NJ/2;
  while(counter<InitNrAg)
    {
      i=(int)(uniform()*10)-5;
      j=(int)(uniform()*10)-5;
      if(i+ii>=NI || i+ii<0 || j+jj>=NJ || j+jj<0)//check whether the selected spot falls within the field
	continue;
      
      if(Ags[ii+i][jj+j]==NULL)
	{
	  if(counter==0)
	    {
	      unfit=TRUE;
	      while(unfit==TRUE)
		{
		  printf("creating agent \n");
		  A=new Agent();
		  A->CreateAgentFromScratch(pagidc_,0);
		  printf("fitness of initial agent %f\n",A->fitness);
		  printf("nonexpfitness of initial agent %f\n",A->nonexpfitness);
		  printf("nrbands of initial agent %i %i\n",A->nrbands,A->nrlongbands);
		  if(A->nrlongbands>=1)//==1)
		    {
		      Ags[ii+i][jj+j]=A;
		      pagidc_++;
		      pnrags_++;
		      unfit=FALSE;
#ifdef RERUN
		      WriteFullAncestry(Ags[ii+i][jj+j]);
#endif		      
		    }
		  else
		    delete A;
		}
	    }
	  else
	    {
	      Arest=new Agent();
	      Arest->CloneAgentFromAgent(A,pagidc_,0);
	      Ags[ii+i][jj+j]=Arest;
	      pagidc_++;
	      pnrags_++;
#ifdef RERUN
	      WriteFullAncestry(Ags[ii+i][jj+j]);
#endif
	    }
	  counter++;
	}
    }
}


void Population::InitPopulationFromGenome(char *fname)
{
  int i,j;
  int ii,jj;
  Agent *A;
  Agent *Arest;
  pagidc_=1;
  // int unfit;
  int counter;
 
  counter=0;
  ii=NI/2;
  jj=NJ/2;
  while(counter<InitNrAg)
    {
      i=(int)(uniform()*10)-5;
      j=(int)(uniform()*10)-5;

      if(Ags[ii+i][jj+j]==NULL)
	{
	  if(counter==0)
	    {
	      //create agent A from genomefile
	      A=new Agent();
	      A->CreateAgentFromGenome(1,1,fname);
	      Ags[ii+i][jj+j]=A;
	      pagidc_++;
	      pnrags_++;
#ifdef RERUN
	      WriteFullAncestry(Ags[ii+i][jj+j]);
#endif		      
	    }
	  else
	    {
	      Arest=new Agent();
	      Arest->CloneAgentFromAgent(A,pagidc_,0);
	      Ags[ii+i][jj+j]=Arest;
	      pagidc_++;
	      pnrags_++;
#ifdef RERUN
	      WriteFullAncestry(Ags[ii+i][jj+j]);
#endif
	    }
	  counter++;
	}
    }
}

void Population::DeathOfAgent(int i,int j)
{
  if(Ags[i][j]!=NULL)
    {
      delete (Ags[i][j]);
      Ags[i][j]=NULL;
      pnrags_--;
    } 
}

void Population::ReproduceAgent(int t,int sourcei,int sourcej,int desti,int destj)
{
  Agent *P;
  Agent *C;
 
  P=Ags[sourcei][sourcej];
  if(P==NULL)
    {
      printf("parent does not exist %i %i\n",sourcei,sourcej);
      exit(1);
    }
  C=new Agent();
  C->CreateAgentFromParent(P,pagidc_,t);
  Ags2[desti][destj]=C;
  pagidc_++;
  pnrags_++;

#ifdef RERUN
  WriteFullAncestry(C);
#endif
}


void Population::Margolus()
{
  int i,j;
  //division type 1
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      Ags2[i][j]=NULL;
  for(i=0;i<NI-1;i+=2)
    for(j=0;j<NJ-1;j+=2)
      {
	if(uniform()<0.5)//clockwise
	  {
	    Ags2[i+1][j]=Ags[i][j];
	    Ags2[i+1][j+1]=Ags[i+1][j];
	    Ags2[i][j+1]=Ags[i+1][j+1];
	    Ags2[i][j]=Ags[i][j+1];
	  }
	else//counterclockwise
	  {
	    Ags2[i][j+1]=Ags[i][j];
	    Ags2[i][j]=Ags[i+1][j];
	    Ags2[i+1][j+1]=Ags[i][j+1];
	    Ags2[i+1][j]=Ags[i+1][j+1];
	  }
      }
  
  //division type 2
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      Ags[i][j]=NULL;
  for(i=0;i<NI;i++)
    {
      j=0;
      Ags[i][j]=Ags2[i][j];
      j=NJ-1;
      Ags[i][j]=Ags2[i][j];
    }
  for(j=0;j<NJ;j++)
    {
      i=0;
      Ags[i][j]=Ags2[i][j];
      i=NI-1;
      Ags[i][j]=Ags2[i][j];
    }
  for(i=1;i<NI-2;i+=2)
    for(j=1;j<NJ-2;j+=2)
      {
	if(uniform()<0.5)//clockwise
	  {
	    Ags[i+1][j]=Ags2[i][j];
	    Ags[i+1][j+1]=Ags2[i+1][j];
	    Ags[i][j+1]=Ags2[i+1][j+1];
	    Ags[i][j]=Ags2[i][j+1];
	  }
	else//counterclockwise
	  {
	    Ags[i][j+1]=Ags2[i][j];
	    Ags[i][j]=Ags2[i+1][j];
	    Ags[i+1][j+1]=Ags2[i][j+1];
	    Ags[i+1][j]=Ags2[i+1][j+1];
	  }
      }
}


void Population::Shuffle()
{
  int i,j;
  int newi,newj;
  int loc;

  //random_shuffle : stl library   
  std::vector<int > locations;
  std::vector<int >::iterator iter;

  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	loc=i*NJ+j;
	locations.push_back(loc);
      }
  std::random_shuffle(locations.begin(),locations.end());

  iter=locations.begin();
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	loc=(*iter);
	newi=loc/NJ;
	newj=loc%NJ;
	Ags2[newi][newj]=Ags[i][j];
	iter++;
      }

  //updating
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      Ags[i][j]=Ags2[i][j];
}



void Population::Step(int t)
{
  int i,j;
  int ii,jj;
  const int RadRep=3;
  static double cumfit[2*RadRep+1][2*RadRep+1];
  double totalfit;
  double fit;
  int memi;
  int memj;
  double temp;
  int nrnb;

  //put Ags2 to zero so no old ones there
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      Ags2[i][j]=NULL;

  //determine who will reproduce 
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]==NULL)
	  {
	    memi=-100;
	    memj=-100;
	   
	    for(ii=0;ii<=2*RadRep;ii++)
	      for(jj=0;jj<=2*RadRep;jj++)
		cumfit[ii][jj]=-1;
	   
	    nrnb=0;
	    totalfit=0;
	    for(ii=max(0,i-RadRep);ii<=min(NI-1,i+RadRep);ii++)
	      for(jj=max(0,j-RadRep);jj<=min(NJ-1,j+RadRep);jj++)
		{		  
		 if(Ags[ii][jj]!=NULL)
		    {
		      fit=Ags[ii][jj]->fitness;
		      if(fit>0)
			{
			  totalfit+=fit;
			  cumfit[ii-i+RadRep][jj-j+RadRep]=totalfit;
			  nrnb++;//count local agents
			}
		    }
		}
	    //check if there are enough agents to reproduce
	    if(nrnb>0 && totalfit>0)
	      {
		//try to select an agent to reproduce
		 temp=uniform()*totalfit;
		 for(ii=0;ii<=2*RadRep;ii++)
		   for(jj=0;jj<=2*RadRep;jj++)
		     {
		       if(temp<cumfit[ii][jj])
			 {
			   memi=ii+i-RadRep;
			   memj=jj+j-RadRep;
			   
			   ii=100;
			   jj=100;
			 }
		     }

		//check if parent selected
		if(memi!=-100)
		  {	    
		    ReproduceAgent(t,memi,memj,i,j);
		  }
	      }
	  }
      }

  //determine who will die
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    if(uniform()<deathrate)
	      DeathOfAgent(i,j);
	    else//if you do not die you go to next step
	      Ags2[i][j]=Ags[i][j];
	  }
      }

  //Update Ags with Ags2
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      Ags[i][j]=Ags2[i][j];

  //!!!!!!!!!!!!!!!!!!
  //Shuffle();
  Margolus();
  //!!!!!!!!!!!!!!!!!! 
}

void Population::WriteFieldToFile(int time)
{
  FILE *f;
  char fname[800],Cname[800];
  int i;
  int j;
  
  //need to make dir in which it all goes
  char dirname[800];
  char command[800];
  
  //printf("am heeeere\n");
  sprintf(dirname,"Field%.10d",time);
  sprintf(command,"mkdir %s/%s",despath,dirname);
  // system(command);
  if(system(command)==-1)//make directory
  { 
    printf("WriteFieldToFile:warning: could not make dir %s. Exiting...\n",dirname);
    exit(1);
  }
  
  
  sprintf(fname,"%s/%s/AgentData",despath,dirname);
  f=fopen(fname,"a"); 
  
  //go through field
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	{
	  fprintf(f,"%i\t%i\t%i\t%f\t%i\n",i,j,Ags[i][j]->agentid,Ags[i][j]->fitness,Ags[i][j]->nrlongbands);
	  Ags[i][j]->WriteGenome(dirname);
	}
    }
 
}

void Population::WritePopFitnessToFile(int time)
{
  FILE *f;
  char fname[800];
  int i;
  int j;
  int nrag;
  double meanfit,maxfit,minfit,stdfit;

  meanfit=0;
  minfit=1000;
  maxfit=-1000;
  stdfit=0;

  nrag=0;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    meanfit+=Ags[i][j]->fitness;
	    stdfit+=Ags[i][j]->fitness*Ags[i][j]->fitness;
	    if(Ags[i][j]->fitness<minfit) 
	      minfit=Ags[i][j]->fitness;
	    if(Ags[i][j]->fitness>maxfit)
	      maxfit=Ags[i][j]->fitness;
	    nrag++;
	  }
    }

  meanfit/=nrag;
  stdfit/=nrag;
  stdfit=stdfit-meanfit*meanfit;
  if(stdfit>0)
    stdfit=sqrt(stdfit);
 
  sprintf(fname,"%s/PopFitnessDynamics",despath);
  f=fopen(fname,"a");
  fprintf(f,"%i\t%i\t%f\t%f\t%f\t%f\n",time,nrag,meanfit,minfit,maxfit,stdfit);
  fclose(f);
}


void Population::WriteGenomeLengthToFile(int time)
{
  FILE *f;
  char fname[400];
  int i;
  int j;
  int nrag;

  double avglen, maxlen,minlen,stdlen,len;
  double avggn,maxgn,mingn,stdgn,gn;
  double avgtfbs, maxtfbs,mintfbs,stdtfbs,tfbs;

  avglen=0;
  avggn=0;
  avgtfbs=0;
  maxlen=0;
  maxgn=0;
  maxtfbs=0;
  minlen=10000;
  mingn=10000;
  mintfbs=10000;
  stdlen=0;
  stdgn=0;
  stdtfbs=0;
  nrag=0;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    len=Ags[i][j]->G->glength_;
	    gn=Ags[i][j]->G->gnrgenes_;
	    tfbs=Ags[i][j]->G->gnrtfbs_;
	    
	    avglen+=len;
	    avggn+=gn;
	    avgtfbs+=tfbs;
	    
	    if(len>maxlen) maxlen=len;
	    if(gn>maxgn) maxgn=gn;
	    if(tfbs>maxtfbs)maxtfbs=tfbs;
	    
	    if(len<minlen) minlen=len;
	    if(gn<mingn) mingn=gn;
	    if(tfbs<mintfbs) mintfbs=tfbs;
	    
	    stdlen+=len*len;
	    stdgn+=gn*gn;
	    stdtfbs+=tfbs*tfbs;

	    nrag++;
	  }
    }

  avglen/=(double)nrag;
  avggn/=(double)nrag;
  avgtfbs/=(double)nrag;

  stdlen/=(double)nrag;
  stdlen=stdlen-avglen*avglen;
  stdlen=sqrt(stdlen);
  stdgn/=(double)nrag;
  stdgn=stdgn-avggn*avggn;
  stdgn=sqrt(stdgn);
  stdtfbs/=(double)nrag;
  stdtfbs=stdtfbs-avgtfbs*avgtfbs;
  stdtfbs=sqrt(stdtfbs);

  sprintf(fname,"%s/PopGenomeLengthDynamics",despath);
  f=fopen(fname,"a");
  fprintf(f,"%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,avggn,mingn,maxgn,stdgn,avgtfbs,mintfbs,maxtfbs,stdtfbs,avglen,minlen,maxlen,stdlen);
  fclose(f);
}


void Population::WriteBandsToFile(int time)
{
  FILE *f;
  char fname[400];
  int i;
  int j;
  int nrag;

  double avgnrlongbands;
  double maxnrlongbands;
  double minnrlongbands;
  double stdnrlongbands;
  double nrlongbands;

  double avgnrbands;
  double maxnrbands;
  double minnrbands;
  double stdnrbands;
  double nrbands;

  avgnrlongbands=0;
  avgnrbands=0;
  maxnrlongbands=0;
  maxnrbands=0;
  minnrlongbands=10000;
  minnrbands=10000;
  stdnrlongbands=0;
  stdnrbands=0;
  // int k;

  nrag=0;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    nrlongbands=Ags[i][j]->nrlongbands;
	    nrbands=Ags[i][j]->nrbands;
	    
            avgnrlongbands+=nrlongbands;
	    avgnrbands+=nrbands;
	   
            if(nrlongbands>maxnrlongbands)
              maxnrlongbands=nrlongbands;
	    if(nrbands>maxnrbands)
	      maxnrbands=nrbands;
	   
            if(nrlongbands<minnrlongbands)
              minnrlongbands=nrlongbands;
	    if(nrbands<minnrbands)
	      minnrbands=nrbands;
	    
            stdnrlongbands+=nrlongbands*nrlongbands;
	    stdnrbands+=nrbands*nrbands;
	  
	    nrag++;
	  }
      }
  avgnrlongbands/=(double)nrag;
  avgnrbands/=(double)nrag;
 
  stdnrlongbands/=(double)nrag;
  stdnrlongbands=stdnrlongbands-avgnrlongbands*avgnrlongbands;
  stdnrlongbands=sqrt(stdnrlongbands);
  stdnrbands/=(double)nrag;
  stdnrbands=stdnrbands-avgnrbands*avgnrbands;
  stdnrbands=sqrt(stdnrbands);
 
  sprintf(fname,"%s/PopBandDynamics",despath);
  f=fopen(fname,"a");
  fprintf(f,"%i\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",time,avgnrbands,minnrbands,maxnrbands,stdnrbands,avgnrlongbands,minnrlongbands,maxnrlongbands,stdnrlongbands);
  fclose(f);
  //1 time
  //2 avgnrbands
  //3 minnrbands
  //4 maxnrbands
  //5 stdnrbands
  //6 avgnrlongbands
  //7 minnrlongbands
  //8 maxnrlongbands
  //9 stdnrlongbands
}


void Population::WriteFittestMatricesToFile(int time)
{
  
  //FILE *f;
  //char fname[600];
  int samplesize;
  int I,J;
  int i,j;
  int selectedags[10000][2];

  for(i=0;i<10000;i++)
    for(j=0;j<2;j++)
      selectedags[i][j]=0;

  //determine fitness of fittest agent
  double maxnonexpfitness=0;
  int Ifittest=-1;
  int Jfittest=-1;
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    if(Ags[i][j]->nonexpfitness>maxnonexpfitness)
	      {
		maxnonexpfitness=Ags[i][j]->nonexpfitness;
		Ifittest=i;
		Jfittest=j;
	      }
	  }
      }
  //store fittest agent
  if(Ifittest>-1 && Jfittest>-1)
    {
      samplesize=0;
      selectedags[samplesize][0]=Ifittest;
      selectedags[samplesize][1]=Jfittest;
      samplesize++;
      //printf("%f %i %i\n",maxnonexpfitness,Ifittest,Jfittest);
    }
  //find other individuals that have the same/similar fitness
  int X,Y,Z;
  int maxsamplesize=10;
  for(i=0;i<1000;i++)
    {
      Z=(int)(uniform()*NI*NJ);
      X=Z/NI;
      Y=Z%NI;

      if(Ags[X][Y]!=NULL)
	{
	  if(Ags[X][Y]->nonexpfitness>=0.95*maxnonexpfitness && samplesize<maxsamplesize)
	    {
	      selectedags[samplesize][0]=X;
	      selectedags[samplesize][1]=Y;
	      samplesize++;
	    }
	}
      if(samplesize==maxsamplesize)
	break;
    }

  //need to make dir in which it all goes
  char dirname[800];
  char command[800];
  sprintf(dirname,"FittestGeneration%.10d",time);
  sprintf(command,"mkdir %s/%s",despath,dirname);
  // system(command);
  if(system(command)==-1)
    { //make directory
      printf("WriteFittestMatricesToFile:warning: could not make dir %s. Exiting...\n",dirname);
      exit(1);
    }

  for(i=0;i<samplesize;i++)
    {
      I=selectedags[i][0];
      J=selectedags[i][1];
      WriteAgentToFile(dirname,I,J);
    }
}

void Population::WriteAncestriesToFile(int t)
{
  FILE *f;
  char fname[800];
  char s[800];
  int i;
  int j;
  iter pi;
  Agent::iter ii;

  //empty population ancestry
  popancestry.clear();

  //copy ancestries of all currently alive and present agents into population ancestry
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    copy(Ags[i][j]->ancestry.begin(),Ags[i][j]->ancestry.end(),std::back_inserter(popancestry));
	    popancestry.push_back(Ags[i][j]->agentid);//agent itself	    
	  }
      }
  //sort population ancestry on number of ancestor and hence temporal order
  popancestry.sort();
  //put non-unique ancestors at end op population ancestry and than remove them
  popancestry.erase(unique(popancestry.begin(),popancestry.end()),popancestry.end()); //why not just use popancestry.unique()?
  
  //write resulting population ancestry with unique ancestors only to file
  sprintf(fname,"%s/PopAncestry%.4d",despath,t);
  f=fopen(fname,"w");
  for(pi=popancestry.begin();pi!=popancestry.end();pi++)
    {
      fprintf(f,"%i ",(*pi));
    }
  fprintf(f,"\n");
  fclose(f);
  sprintf(s,"gzip %s",fname);
  if(system(s)==-1)
    { //zip data
      printf("WriteAncestriesToFile:warning: could not zip %s. Exiting...\n",fname);
      exit(1);
    }
  
  //write a list of ancestries of all currently present agents to file
  sprintf(fname,"%s/Ancestry%.4d",despath,t);
  f=fopen(fname,"w");
  for(i=0;i<NI;i++)
    for(j=0;j<NJ;j++)
      {
	if(Ags[i][j]!=NULL)
	  {
	    for(ii=(*Ags[i][j]).ancestry.begin();ii!=(*Ags[i][j]).ancestry.end();ii++)
	      {
		fprintf(f,"%i ",(*ii));
	      }
	    fprintf(f,"%i ",Ags[i][j]->agentid);
	    fprintf(f,"\n");
	  }
      }

  fclose(f);
  sprintf(s,"gzip %s",fname);
  // system(s);
  if(system(s)==-1)
    { //zip data
      printf("WriteAncestriesToFile:warning 2: could not zip %s. Exiting...\n",fname);
      exit(1);
    }
  

}


void Population::ReadAncestorsFromFile(char *fname)
{
  FILE *f;
  char fulfname[800];
  int intdummy;

  //sprintf(fulfname,"%s/%s",despath,fname);
  f=fopen(fname,"r");
  while(fscanf(f,"%i",&intdummy)!=EOF)
    {
      popancestry.push_back(intdummy);
    }
  fclose(f);
}


void Population::WriteFullAncestry(Agent *C)//RERUN
{
  FILE *f;
  char fname[800];

  //printf("agent id %i\n",C->agentid);


  if(!popancestry.empty())
    {
      if(find(popancestry.begin(),popancestry.end(),C->agentid)!=popancestry.end())//agent is important for ancestry
	{
	  //printf("saved \n");
	  sprintf(fname,"%s/FullAncestry",despath);
	  f=fopen(fname,"a");
	  fprintf(f,"%i %i %i %i %i %f %i %i %i\n",C->agentid,C->tbirth,C->nrbands,C->nrlongbands,C->cells.size(),C->fitness,C->G->glength_,C->G->gnrgenes_,C->G->gnrtfbs_);
	  fclose(f);

	  // sprintf(fname,"full");
	  C->WriteGenome("AncestryGenomes");

	  popancestry.remove(pagidc_);
	}
    }
}



void Population::WriteAgentToFile(char *dirname,int i,int j)
{
  Ags[i][j]->WriteEmbryology(dirname);
  Ags[i][j]->WriteDivisionProfile(dirname);
  Ags[i][j]->WriteGenome(dirname);
  //Ags[i][j]->WriteSignalProfiles(dirname);
  Ags[i][j]->WriteFitnessDetails(dirname);
}

