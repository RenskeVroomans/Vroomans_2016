#include "Header.hh"
#include "Agent.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Network.hh"
#include "Edge.hh"
#include "Vertex.hh"
#include "Misc.hh"
#include "dSFMT.h"
#include <sys/stat.h>

dsfmt_t dsfmt;
double gauss[10000];

char writepath[500];
char filename[500];

int FullAnalysis;

void WriteVertexList(Agent *A,int I)
{
  FILE *f1;
  FILE *f2;
  char s1[800];
  char s2[800];
  int i,j;
  int c1;
  int c2;
  int totalcount;
  Network::iterv iv;
  Network::iterel el;
  Network::itere e;
  int genelist[100];
  int genesused[NrGeneTypes+1];
  int genetype;
  int Adj[100][100];
  Gene *gt;
  Gene *gt2;

  for(i=0;i<100;i++)
    genelist[i]=-1;
  for(i=0;i<NrGeneTypes+1;i++)
    genesused[i]=0;
  for(i=0;i<100;i++)
    for(j=0;j<100;j++)
      Adj[i][j]=0;

  //first make a list of the different genes in the genome and their ordering
  //make sure to distinghuish multiple copies of the same gene type as different genes
  i=0;//count vertex number = nr of gene (not gene type nr)
  iv=A->N->VL->begin();
  while(iv!=A->N->VL->end())
    {
      gt=(*iv)->Gen;
      //genetype is increased by 1 so 0->1 n->n+1
      genetype=gt->type+1;
      //first time gene of certain type is encountered
      if(genesused[genetype]==0)
	{
	  genelist[i]=genetype;//store that i-th gene in list is of this type
	  genesused[genetype]=1;//store that a gene of this type has been encountered
	}
      else
	{
	  //if not first copy of gene: i.e genesused is no longer 0
	  //make sure that you use a different number for new copies of the same gene type
	  //such that you know that it is a different gene /node but of the same type
	  //second copy: genetypenr*100
	  //third copy: genetypenr*100*100 
	  genelist[i]=genetype*pow(100.,(double)genesused[genetype]);
	  genesused[genetype]++;
	}
      i++;
      iv++;
    }
  printf("nr of genes in the genome %i\n",i);
  int nrgenes=i;

  if(I==1)
    sprintf(s2,"%s/GeneOrderListAgent%.10doriginal",despath,A->agentid);
  else if(I==2)
    sprintf(s2,"%s/GeneOrderListAgent%.10dcore",despath,A->agentid);
  f2=fopen(s2,"w");
  for(i=0;i<nrgenes;i++)
    fprintf(f2,"%i\t%i\n",i,genelist[i]);
  fclose(f2);
   
  
  if(I==1)
    sprintf(s2,"%s/VertexListAgent%.10doriginal",despath,A->agentid);
  else if(I==2)
    sprintf(s2,"%s/VertexListAgent%.10dcore",despath,A->agentid);
  f2=fopen(s2,"w");

  //now go through the adjacency list to determine the regulatory links between genes
  //c1 and c2 count the number/position the gene has in the vertex list/genome
  //c1 monotonically walks along the to be regulated genes in the adjacency/vertex list
  //c2 is the position of a certain gene that has an incoming regulatory link to gene c1
  //in the vertex list, so now we know there is a link from gene gt2 to gene gt and hence
  //from gene number c2 to gene number c1
  el=A->N->AL->begin();
  c1=0;
  totalcount=0;
  while(el!=A->N->AL->end())
    {
      for(e=el->begin();e!=el->end();e++)
	{
	  gt=(*e)->V->Gen;//gene that we found in adjacency list: incoming regulatory link
	  c2=0;
	  iv=A->N->VL->begin();
	  while(iv!=A->N->VL->end())
	    {
	      gt2=(*iv)->Gen;//gene that we found in vertex list: genes to be regulated
	      //is it correct that this compares the gene-objects,
	      //thus finding differences between different genes of 
	      //the same type?
	      //hmm, I guess it compares the pointers to the genes
	      //and different genes of the same type would still
	      //have different pointers to them, so yes
	      if(gt==gt2)//are the genes from adjacency and vertex list the same
		{        
		  Adj[c1][c2]++;
		  Adj[c2][c1]++;
		  fprintf(f2,"%i %i\n",genelist[c2],genelist[c1]);
		  totalcount++;
		}
	      iv++;
	      c2++;
	    }
	}
      c1++;
      iv++;
      el++;
    }
  fclose(f2);
  printf("nr of connections in the network %i\n",totalcount);


  if(I==1)
    {
      int nrinputs=0;
      iv=A->N->VL->begin();
      el=A->N->AL->begin();
      while(iv!=A->N->VL->end())
	{
	  if((*iv)->Gen->type==SegmGeneNr)
	    {
	       for(e=el->begin();e!=el->end();e++)
		 {
		   nrinputs++;
		 }
	    }
	  iv++;
	  el++;
	}
      
      sprintf(s1,"%s/InputsToSegmentationGene",despath);
      f1=fopen(s1,"a");
      fprintf(f1,"%i\t",A->agentid);
      fprintf(f1,"%i\n",nrinputs);
      fclose(f1);
    }
}

Agent* PruneAgent(Agent *A,int c)
{
  Agent *Ad;
  Agent *Ap;
  Agent *Ap2;
  int PruneIter=200;

  printf("%i\n",c);

  if(c== 2 || c==3|| c==5 || c==7 || c==9)
  {
    printf("invalid option for PruneAgent: not checking domains now. exiting...\n");
    exit(1);
  }
  
  if(!(c==7 || c==10))//if(c!=7)
    {
      Ap=new Agent();
      Ap->CloneAgentFromAgent(A,A->agentid,0);
      
      int counter=0;
      while(counter<PruneIter)
	{
	  Ad=new Agent();
	  Ad->InitAgent(A->agentid,0,0);
	  Ad->G->CloneGenome(Ap->G);
	  Ad->G->PruneGenome();
	  Ad->N->BuildNetwork(Ad->G);
	  Ad->DevelopAgent();
	  Ad->DetermineFitness();
	  if(Ad->CompareAgents(A,c)==0)//ok
	    {
	      delete Ap;
	      Ap=Ad;
	      Ad=NULL;
	    }
	  else
	    {
	      delete Ad;
	      Ad=NULL;
	    }
	  counter++;
	}
    }
  else
    {
      //in deze komt t resultaat
      Ap=new Agent();
      Ap->InitAgent(A->agentid,0,0);
      Ap->G->CloneGenome(A->G);
      Ap->G->RemoveSegmentationGene();
      Ap->N->BuildNetwork(Ap->G);
      Ap->DevelopAgent();
      Ap->DetermineFitness();

      if(c==7)
	{
      
	  //deze gebruik je om steeds mee te vergelijken
	  Ap2=new Agent();
	  Ap2->CloneAgentFromAgentForSafety(Ap,Ap->agentid,0);
	  
	  int counter=0;
	  while(counter<PruneIter)
	    {
	      Ad=new Agent();
	      Ad->InitAgent(A->agentid,0,0);
	      Ad->G->CloneGenome(Ap->G);
	      Ad->G->PruneGenome();
	      Ad->N->BuildNetwork(Ad->G);
	      Ad->DevelopAgent();
	      Ad->DetermineFitness();
	      if(Ad->CompareAgents(Ap2,c)==0)//ok
		{
		  delete Ap;
		  Ap=Ad;
		  Ad=NULL;
		}
	      else
		{
		  delete Ad;
		  Ad=NULL;
		}
	      counter++;
	    }
	  
	  delete Ap2;
	}
    }

  //Ap->WriteEmbryology(c);
  Ap->WriteGenome(c);
  Ap->WriteNetwork(c);
  Ap->WriteExtraSignalProfiles(c);

  printf("!!!\n");
  return Ap;
}

void Start(int argc,char **argv)
{
  
  struct stat sb;        //allows to check for the existence of directories
  
  if(argc<5)
  {
    printf("usage: <program name> <readpath> <destination path> <inputfile> <FULL ->(0/1)>\n");
    exit(1);
  }
  else
  {
    strcpy(despath,argv[1]);
    strcpy(writepath,argv[2]);
    strcpy(filename,argv[3]);
    FullAnalysis=atoi(argv[4]);
  }
  
  char makedir[]="mkdir ";
  char makecommand[50];
  strcpy(makecommand, makedir);
  strcat(makecommand,writepath);
  
  /// make directory if necessary ///
  if (stat(writepath, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists
      //do nothing;
  }
  
  else if(stat(writepath, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
      printf("given destination is not a directory. give another name\n");
      exit(1);
  }
  else{ //directory does not exist: make it.
      if(system(makecommand)==-1){ //make directory for data
	printf("warning: could not make directory %s. Exiting...\n",writepath);
	exit(1);
      }
  }
  
  
  
}

int main(int argc, char **argv)
{
  FILE *f;
  FILE *f2;
  FILE *f3;
  FILE *f4;
  char fname[1000];
  char fname2[1000];
  char fname3[1000];
  char fname4[1000];
  Agent *A1;
  Agent *A2;
  int temp;
  int id1;
  int id2;
  int genosim;
  double genodist;
  double phenodist;
  double cumgenodist=0;
  double cumphenodist=0;
  Agent *Acore;
  Agent *Acoreminsegm;
  Agent *Anrsegm;
  Agent *Asegm;
  Agent *Aminsegm;
  
  /** initialize random numbers */
  fillgauss();
  dsfmt_init_gen_rand(&dsfmt,94);//dsfmt_init_gen_rand(&dsfmt,11);
  
  /**read command line parameters */
  Start(argc,argv);
  
  /**open the input file: contains just the agent IDs: ancestor file */
  sprintf(fname,"%s",filename);
  f=fopen(fname,"r");
  /**open the output file into the same directory as the input file */
  sprintf(fname2,"%s/%s%s",writepath,filename,"compressed");
  f2=fopen(fname2,"a");
  
  /** read (and print) ID of first agent */
  fscanf(f,"%i",&temp);//first nr =-1, no agent
  fscanf(f,"%i",&id1);
  fprintf(f2,"%i ",temp);//first nr =-1, no agent
  fprintf(f2,"%i ",id1);
  fclose(f2);
  
  /************************************/
  /** start analysing the first agent */
  A1=new Agent();
  
  A1->CreateAgentFromFile(despath,id1); //develop it anew
  
  //A1->WriteExtensiveEmbryology(0); //make more pictures
  //A1->WriteGenome(0);
  A1->WriteNetwork(0);    
  //A1->WriteSignalProfiles(0);
  //single file for ancestry
  A1->WriteBasicProperties(A1->agentid);
  A1->DetermineGenomeAndNetworkProperties();
  A1->WriteGenomeAndNetworkProperties(A1->agentid);
  A1->DetermineAttractorProperties();
  A1->WriteAttractorProperties(A1->agentid);
  A1->DetermineLoopAndMotifProperties();
  A1->WriteLoopAndMotifProperties(A1->agentid);
  
  if(FullAnalysis){
  sprintf(fname3,"%s/%s",writepath,"GeneTFBSConnectionNumbers");
  f3=fopen(fname3,"a");
  fprintf(f3,"%i\t",A1->agentid);
  fprintf(f3,"%i\t%i\t",A1->G->gnrgenes_,A1->G->gnrtfbs_);// A1->N->nrvertices

  Acore=PruneAgent(A1,1);
  fprintf(f3,"%i\t%i\t",Acore->G->gnrgenes_,Acore->G->gnrtfbs_);//5,6,7 
  
  Asegm=PruneAgent(Acore,8);
  fprintf(f3,"%i\t%i\t",Asegm->G->gnrgenes_,Asegm->G->gnrtfbs_);//8,9,10

  Anrsegm=PruneAgent(Asegm,4);
  fprintf(f3,"%i\t%i\t",Anrsegm->G->gnrgenes_,Anrsegm->G->gnrtfbs_);//14,15,16

  Aminsegm=PruneAgent(Asegm,6);
  fprintf(f3,"%i\t%i\t",Aminsegm->G->gnrgenes_,Aminsegm->G->gnrtfbs_);//20,21,22

  Acoreminsegm=PruneAgent(Acore,10);
  fprintf(f3,"%i\t%i\t",Acoreminsegm->G->gnrgenes_,Acoreminsegm->G->gnrtfbs_);//32,33,34

  fprintf(f3,"\n");

  fclose(f3);
  
  WriteVertexList(A1,1);
  WriteVertexList(Acore,2);
  
  sprintf(fname4,"%s/SegmentNumbers",writepath);
  f4=fopen(fname4,"a");
  fprintf(f4,"%i\t",A1->agentid);//1
  fprintf(f4,"%i\t%i\t",A1->nrlongbands,A1->nrlonganddifbands);//2,3
  fprintf(f4,"%i\t%i\t",Acore->nrlongbands,Acore->nrlonganddifbands);//4,5
  fprintf(f4,"%i\t%i\t",Acoreminsegm->nrlongbands,Acoreminsegm->nrlonganddifbands);//6,7
  fprintf(f4,"%i\t%i\t",Aminsegm->nrlongbands,Aminsegm->nrlonganddifbands);//8,9
  fprintf(f4,"%i\t%i\t",Acore->nrlongbands-Aminsegm->nrlongbands,Acore->nrlonganddifbands-Aminsegm->nrlonganddifbands);
  fprintf(f4,"%f\t%f\t",(double)Aminsegm->nrlongbands/(double)Acore->nrlongbands,(double)Aminsegm->nrlonganddifbands/(double)Acore->nrlonganddifbands);
  fprintf(f4,"\n");
  fclose(f4);
  
  delete Acore;
  delete Anrsegm;
  delete Aminsegm;
  delete Asegm;
  delete Acoreminsegm;
  } //fullanalysis
  else
    WriteVertexList(A1,1);
  
  
  /*********************************/
  /** now move to the other agents */
  
   while(fscanf(f,"%i",&id2)!=EOF)
    {
      A2=new Agent();
      //only initialise the agent from file
      A2->InitAgentFromFile(despath,id2);
      //determine whether it has same or different genome as previous agent
      A1->DetermineGenotypeSimilarity(A2,&genosim);
      if(genosim==0)//genome is same as that of previous agent
	{
	  //no need to compute development of agent A2
	  //as it will be exactly the same as those of agent A1
	  //use properties of agent A1 to describe agent A2
	  //single agent files: embryology, genome, network 
	  //signal profiles: no need to do again
	  //whole ancestry files: use agent A1 properties to
	  //write agent A2 properties
	  A1->WriteBasicProperties(A2->agentid);
	  A1->WriteGenomeAndNetworkProperties(A2->agentid);
	  //A1->WriteAttractorProperties(A2->agentid);
	  //A1->WriteLoopAndMotifProperties(A2->agentid);
	  A1->WriteGenoNetworkPhenoDistance(A2->agentid,genodist,phenodist,cumgenodist,cumphenodist);
	  delete (A2);
	  A2=NULL; 
	}
      else//genome is different from that of previous agent
	{ 
	  //so we need to fully develop this agent to determine properties
	  A2->FinishAgentFromFile();
	  //single file per agent
	  A2->WriteEmbryology(0);
	  //A2->WriteGenome(0);
	  A2->WriteNetwork(0);
	  //A2->WriteSignalProfiles(0);
	  //single file for ancestry
	  A2->WriteBasicProperties(A2->agentid);
	  A2->DetermineGenomeAndNetworkProperties();
	  A2->WriteGenomeAndNetworkProperties(A2->agentid);
	  A2->DetermineAttractorProperties();
	  A2->WriteAttractorProperties(A2->agentid);
	  A2->DetermineLoopAndMotifProperties();
	  A2->WriteLoopAndMotifProperties(A2->agentid);
	  A2->DetermineGenoNetworkPhenoDistance(A1,&genodist,&phenodist);
	  cumgenodist+=genodist;
	  cumphenodist+=phenodist;
	  A2->WriteGenoNetworkPhenoDistance(A2->agentid,genodist,phenodist,cumgenodist,cumphenodist);
	  

	  #ifdef HOMOGENEOUSSEGMENTS
	  if(A1->nrhomogbands!=A2->nrhomogbands || A1->nrlonganddifbands!=A2->nrlonganddifbands)
	    #else
	  if(A1->nrlongbands!=A2->nrlongbands || A1->nrlonganddifbands!=A2->nrlonganddifbands)
	    #endif  
	  {
	    f2=fopen(fname2,"a");
	    fprintf(f2,"%i ",A2->agentid);
	    //printf("%i\n",A2->agentid);
	    fclose(f2);
	  }
	  /*******************************************************************/
	  
	  if(FullAnalysis){
	  f3=fopen(fname3,"a");
	  
	  fprintf(f3,"%i\t",A2->agentid);
	  fprintf(f3,"%i\t%i\t",A2->G->gnrgenes_,A2->G->gnrtfbs_);//2,3,4
	  
	  Acore=PruneAgent(A2,1);
	  fprintf(f3,"%i\t%i\t",Acore->G->gnrgenes_,Acore->G->gnrtfbs_);//5,6,7
	  
	   
	  Asegm=PruneAgent(Acore,8);
	  fprintf(f3,"%i\t%i\t",Asegm->G->gnrgenes_,Asegm->G->gnrtfbs_);//8,9,10
	  
	  Anrsegm=PruneAgent(Asegm,4);
	  fprintf(f3,"%i\t%i\t",Anrsegm->G->gnrgenes_,Anrsegm->G->gnrtfbs_);//14,15,16
	  
	  Aminsegm=PruneAgent(Asegm,6);
	  fprintf(f3,"%i\t%i\t",Aminsegm->G->gnrgenes_,Aminsegm->G->gnrtfbs_);//20,21,22
	  
	  Acoreminsegm=PruneAgent(Acore,10);
	  fprintf(f3,"%i\t%i\t",Acoreminsegm->G->gnrgenes_,Acoreminsegm->G->gnrtfbs_);//32,33,34
	  
	  fprintf(f3,"\n");

	  fclose(f3);

	  WriteVertexList(A2,1);
	  WriteVertexList(Acore,2);


	  sprintf(fname4,"%s/SegmentAndDomainNumbers",writepath);
	  f4=fopen(fname4,"a");
	  fprintf(f4,"%i\t",A2->agentid);
	  fprintf(f4,"%i\t%i\t",A2->nrlongbands,A2->nrlonganddifbands);//2,3
	  fprintf(f4,"%i\t%i\t",Acore->nrlongbands,Acore->nrlonganddifbands);//4,5
	  fprintf(f4,"%i\t%i\t",Acoreminsegm->nrlongbands,Acoreminsegm->nrlonganddifbands);//6,7
	  fprintf(f4,"%i\t%i\t",Aminsegm->nrlongbands,Aminsegm->nrlonganddifbands);//8,9
	  fprintf(f4,"%i\t%i\t",Acore->nrlongbands-Aminsegm->nrlongbands,Acore->nrlonganddifbands-Aminsegm->nrlonganddifbands);//12,13
	  fprintf(f4,"%f\t%f\t",(double)Aminsegm->nrlongbands/(double)Acore->nrlongbands,(double)Aminsegm->nrlonganddifbands/(double)Acore->nrlonganddifbands);//15,17
	  fprintf(f4,"\n");
	  fclose(f4);


	  //printf("overlap %i %i unique %i %i frac %f %f\n",nroverlapgenes,nroverlapconnections,nruniquegenes,nruniqueconnections,(double)nroverlapgenes/(double)nruniquegenes,(double)nroverlapconnections/(double)nruniqueconnections);
	  
	  delete Acore;
	  delete Anrsegm;
	  delete Aminsegm;
	  delete Asegm;
	  delete Acoreminsegm;
	  
	  }//fullanalysis
	  
	  else
	    WriteVertexList(A2,1);
	  
	  /*******************************************************************/
	  

	  delete (A1);
	  A1=NULL;
	  A1=A2;
	  A2=NULL;
	}
    }
  f2=fopen(fname2,"a");
  fprintf(f2,"%i ",A1->agentid);
  fprintf(f2,"\n");
  fclose(f2);
  //printf("%i\n",A1->agentid);
  fclose(f);
}