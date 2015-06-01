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
#include "IO.hh"
#include <sys/stat.h>
#include <png.h>

dsfmt_t dsfmt;
double gauss[10000];

int AgentID;
int extravid;
char filename[500];
char writepath[500];
char readpath[500];
char extrapath[500];//for networks
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

  int accepted=0, nonaccepted=0;
  
  printf("%i\n",c);

  if(c== 2 || c==3|| c==5 || c==7 || c==9 ||c==10 )
  {
    printf("invalid option for PruneAgent: not checking domains now. exiting...\n");
    exit(1);
  }
  
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
    Ad->DetermineFitness(1);
    //printf("new agent: %d bands, original: %d\n", Ad->nrlongbands, A->nrlongbands);
    
    if(Ad->CompareAgents(A,c)==0)
    {
      delete Ap;
      Ap=Ad;
      Ad=NULL;
     // printf("prune accepted.\n\n");
      accepted++;
    }
    else
    {
      delete Ad;
      Ad=NULL;
      //printf("prune not accepted.\n\n");
      nonaccepted++;
    }
    counter++;
  }
  
  //Ap->WriteEmbryology(c);
  Ap->WriteGenome(c);
  Ap->WriteNetwork(c);
  Ap->WriteExtraSignalProfiles(c);
  Ap->WriteSegmentEmbryology(c);
  
  printf("accepted:%d, nonaccepted: %d\n",accepted, nonaccepted);
  printf("!!!\n");
  return Ap;
}

void WriteStripePatterns(int **array,int nrnetworks,int dsize)
{
  int i,j,k,l;
  int ii,jj;
  FILE *PNGFileP;  
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep row_pointer;
  char fname[500];
  const int zoom=8;
  const int WW=zoom*dsize;
  const int LL=zoom*(nrnetworks+nrnetworks-1);
  int celltypes[LL][WW];
  unsigned char RGBdata[LL*WW*3];
  int type;
  int color;

  for(i=0;i<nrnetworks;i++)
  {
    for(j=0;j<dsize;j++)
    {
      if(array[i][j]>ThOn)
	color=0;
      else
	color=1;
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<zoom;jj++)
	  celltypes[(zoom)*(i*2)+ii][zoom*j+jj]=color;
    }
    if(i<nrnetworks-1){
      for(ii=0;ii<zoom;ii++)
	for(jj=0;jj<WW;jj++)
	  celltypes[zoom*(i*2+1)+ii][jj]=1;
    }
  }

    //black/ white scale
    for(i=0;i<LL;i++){
      for(j=0;j<WW;j++)
      {
	RGBdata[i*WW*3+j*3+0]=celltypes[i][j]*255;
	RGBdata[i*WW*3+j*3+1]=255*celltypes[i][j];
	RGBdata[i*WW*3+j*3+2]=255*celltypes[i][j];
      }
    }
    printf("so far in pictures\n");
    
    sprintf(fname,"%s/Segments%.10d_allvariants.png",despath,AgentID);
    
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
  // write header
  png_write_info ( png_ptr, info_ptr );
  // write out image, one row at a time 
  int row;
  for ( row = 0; row <LL; row++) 
    {
      row_pointer = ( RGBdata + WW * row * 3 );
      png_write_rows ( png_ptr, &row_pointer, 1 );
    }
  // flush all info to file 
  png_write_end ( png_ptr, info_ptr );
  fflush ( PNGFileP );
  png_destroy_write_struct ( &png_ptr,&info_ptr);
  fclose(PNGFileP);

  //printf("spacetimeplot development png written to file\n");

}

static char** duplicateArgv(int argc, char **argv, int offset)                                                                                                                                                         
{ 
  char **array;
  
  //const char name[]="programname";
  
  array = (char **)malloc(sizeof(char*)*(argc-offset+2)); // one for program name plus one extra which will mark the end of the array
  int index;
  
  array[0] = (char *)malloc( strlen(argv[0])+1 ); // add one for the \0
  strcpy(array[0],argv[0]);	  
  
  for(index=0; index<offset; index++)
    ++argv;
  
  for (index = 1; index < argc-offset+1; ++index)
  {
    /*if(index==0){
      array[index] = (char *)malloc( strlen(argv[0])+1 ); // add one for the \0
      strcpy(array[index],argv[0]);
    }*/
    //else{
      array[index] = (char *)malloc( strlen(*argv)+1 ); // add one for the \0
      strcpy(array[index], *argv);
      ++argv;
    //}
  }
  array[index] = NULL; // end of array so later you can do while (array[i++]!=NULL)
  
  
  return array;
}

void Start(int argc,char **argv)
{
  
  struct stat sb;        //allows to check for the existence of directories
  char **argv2;
  
  if(argc<6)
  {
    printf("usage: <program name> <agent number> <writepath> <extra Cell network video> <FULL ->(0/1)> <generic options (see --help)\n");
    exit(1);
  }
  else
  {
    AgentID=atoi(argv[1]);
    strcpy(writepath,argv[2]);
    extravid=atoi(argv[3]);
    FullAnalysis=atoi(argv[4]);
  }
  
  argv2=duplicateArgv(argc,argv,5);
  
  //  printf("full analysis: %d\n",FullAnalysis);
  ReadPars(argc-4,argv2);// pass a selection of the command line to readpars
  
  strcpy(readpath,despath); //the path that you give as despath now, is the path to read from (saves a program option)
  strcpy(despath, writepath); //now, we need to write to despath (elsewhere in one of the functions)
  
  char makedir[]="mkdir ";
  char makecommand[100];
  strcpy(makecommand, makedir);
  strcat(makecommand,despath);
  
 
  

  char testagent[100];
  strcpy(testagent, readpath);
  strcat(testagent, "CodedGenomeAgent000");
  strcat(testagent, argv[1]);
  if(stat(testagent, &sb) != 0){
    printf("agent does not exist in this run. Exiting...\n");
    exit(1);
  }
  
  
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
  
  ///make directory for movies if necessary ///
  strcpy(extrapath,writepath);
  strcat(extrapath, "/networkmovies");
  strcat(makecommand,"/networkmovies");
  
  if (stat(extrapath, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists
      //do nothing;
  }
  
  else if(stat(extrapath, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
      printf("destination for movies is already a file. give another name\n");
      exit(1);
  }
  else{ //directory does not exist: make it.
      if(system(makecommand)==-1){ //make directory for data
	printf("warning: could not make directory %s. Exiting...\n",extrapath);
	exit(1);
      }
  }
 
    
}

int main(int argc, char **argv)
{
  FILE *f3;
  FILE *f4;
 
  char fname3[1000];
  char fname4[1000];
  Agent *A1;
  Agent *Acore;
  Agent *Aminsegm;
  Agent *Anrsegm;
  Agent *Asegm;
    
    
  /**read command line parameters */
  Start(argc,argv);
 
  /** initialize random numbers */
  fillgauss();
  dsfmt_init_gen_rand(&dsfmt,seedinitpop);//dsfmt_init_gen_rand(&dsfmt,11);
  
  /************************************/
  /** start analysing the first agent */
  A1=new Agent();

  A1->CreateAgentFromFile(readpath,AgentID); //develop it anew
  A1->WriteSegmentEmbryology(0); //make more pictures
  A1->WriteEmbryology("");
  A1->WriteNetwork(0);    
  A1->WriteExtraSignalProfiles(0);
  //single file for ancestry
  A1->WriteBasicProperties(A1->agentid);
  A1->DetermineGenomeAndNetworkProperties();
  A1->WriteGenomeAndNetworkProperties(A1->agentid);
   A1->DetermineAttractorProperties();
  A1->WriteAttractorProperties(A1->agentid);
  A1->DetermineLoopAndMotifProperties();
  A1->WriteLoopAndMotifProperties(A1->agentid);
  
  A1->WriteNetworkProgression(0,extrapath,0);
  A1->WriteNetworkProgression(0,extrapath,10);
  if (extravid)
    A1->WriteNetworkProgression(0,extrapath,extravid);
  
  if(FullAnalysis){
    
    int *oscarray;//[NrGeneTypes];
    oscarray=(int *)calloc((size_t)NrGeneTypes, sizeof(int ));
    if(oscarray==NULL)
    {
      printf("error in memory allocation segarray. Abort...\n");
      exit(1);
    }
    
    sprintf(fname4,"%s/%s",writepath,"OscillatingProteins");
    f4=fopen(fname4,"a");
  
    sprintf(fname3,"%s/%s",writepath,"GeneTFBSConnectionNumbers");
    f3=fopen(fname3,"a");
    
    /** original **/
    fprintf(f3,"%i\t",A1->agentid);
    fprintf(f3,"%i\t%i\t",A1->G->gnrgenes_,A1->G->gnrtfbs_);// A1->N->nrvertices
    
    A1->FindProteinOscillations(oscarray);
    
    fprintf(f4,"original:\t");
    for(int i=0;i<NrGeneTypes; i++)
    {
      if(oscarray[i])
	fprintf(f4,"%i\t",i);
    }
    fprintf(f4,"\n");
    
    /**core network **/
    Acore=PruneAgent(A1,1);
    fprintf(f3,"%i\t%i\t",Acore->G->gnrgenes_,Acore->G->gnrtfbs_);//5,6,7 
    Acore->FindProteinOscillations(oscarray);
    
    fprintf(f4,"core:\t");
    for(int i=0;i<NrGeneTypes; i++)
    {
      if(oscarray[i])
	fprintf(f4,"%i\t",i);
    }
    fprintf(f4,"\n");
    
     /**segmentation network **/
     Asegm=PruneAgent(A1,8);//PruneAgent(Acore,8);
     fprintf(f3,"%i\t%i\t",Asegm->G->gnrgenes_,Asegm->G->gnrtfbs_);//8,9,10
     Asegm->FindProteinOscillations(oscarray);
     fprintf(f4,"segm:\t");
     for(int i=0;i<NrGeneTypes; i++)
     {
       if(oscarray[i])
	 fprintf(f4,"%i\t",i);
     }
     fprintf(f4,"\n");
     Asegm->WriteNetworkProgression(8,extrapath,0);
     Asegm->WriteNetworkProgression(8,extrapath,10);
     if (extravid)
       Asegm->WriteNetworkProgression(8,extrapath,extravid);
     /**nr segments network **/
    Anrsegm=PruneAgent(A1,4);//////PruneAgent(Asegm,4);
    fprintf(f3,"%i\t%i\t",Anrsegm->G->gnrgenes_,Anrsegm->G->gnrtfbs_);//14,15,16
    Anrsegm->FindProteinOscillations(oscarray);
    fprintf(f4,"nrsegm:\t");
    for(int i=0;i<NrGeneTypes; i++)
    {
      if(oscarray[i])
	fprintf(f4,"%i\t",i);
    }
    fprintf(f4,"\n");
    
    Anrsegm->WriteNetworkProgression(4,extrapath,0);
    Anrsegm->WriteNetworkProgression(4,extrapath,10);
    if(extravid)
      Anrsegm->WriteNetworkProgression(4,extrapath,extravid);
    
    /**min segments network **/
    Aminsegm=PruneAgent(Anrsegm,6);////PruneAgent(A1,6);
    fprintf(f3,"%i\t%i\t",Aminsegm->G->gnrgenes_,Aminsegm->G->gnrtfbs_);//20,21,22
    Aminsegm->FindProteinOscillations(oscarray);
    fprintf(f4,"minsegm:\t");
    for(int i=0;i<NrGeneTypes; i++)
    {
      if(oscarray[i])
	fprintf(f4,"%i\t",i);
    }
    fprintf(f4,"\n");
    
    fprintf(f3,"\n");

    Aminsegm->WriteNetworkProgression(6,extrapath,0);
    Aminsegm->WriteNetworkProgression(6,extrapath,10);
    if (extravid)
      Aminsegm->WriteNetworkProgression(6,extrapath,extravid);
    fclose(f3);
    fclose(f4);
    //WriteVertexList(A1,1);
    //WriteVertexList(Acore,2);
    
    //create combinatorial picture with segments of all reduced agents
    int **segarray;
    int networks=5;
   
    segarray=(int **)calloc((size_t)networks, sizeof(int *));
    if(segarray==NULL)
    {
      printf("error in memory allocation segarray. Abort...\n");
      exit(1);
    }
    segarray[0]=(int *)calloc((size_t)networks*NrFinalCells, sizeof(int));
    if(segarray[0]==NULL)
    {
      printf("error in memory allocation segarray. Abort...\n");
      exit(1);
    }
    for(int i=1;i<networks;i++)
      segarray[i]=segarray[i-1]+NrFinalCells;
   
    int capture=(int)(NrDevSteps/StorageInt);
   
    for(int i=0; i<NrFinalCells; i++)
      segarray[0][i]=A1->E[capture][i][5];
    for(int i=0; i<NrFinalCells; i++)
      segarray[1][i]=Acore->E[capture][i][5];
    for(int i=0; i<NrFinalCells; i++)
      segarray[2][i]=Asegm->E[capture][i][5];
    for(int i=0; i<NrFinalCells; i++)
      segarray[3][i]=Anrsegm->E[capture][i][5];
    for(int i=0; i<NrFinalCells; i++)
      segarray[4][i]=Aminsegm->E[capture][i][5];
    
    
    printf("so far so good!\n");
    WriteStripePatterns(segarray, networks,NrFinalCells);
    
    
    delete Acore;
    delete Anrsegm;
    delete Aminsegm;
    delete Asegm;
    
   
  } //fullanalysis
  
  else
    WriteVertexList(A1,1);
  
  
  
  return 0;
  
}
