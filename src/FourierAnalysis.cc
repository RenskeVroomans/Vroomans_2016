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
#include <fftw3.h>
#include <deque>

#define FIXCELLS 120 //nr of cells for the analysis with fixed concentration of morphogen
#define DECCELLS 50 //nr of cells for the analysis with decaying concentration of morphogen
#define EXTRATIME 600 //extra time for the fixed-concentration cells to run (to do joint time freq analysis)
#define WINDOWSIZE 300 //size of window for sliding-window analysis
dsfmt_t dsfmt;
double gauss[10000];

char genomefile[500];
char writepath[500];
char extrapath[500];//for networks

vector<vector<vector <double> > > store; //genes, freq series, real/im part


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
    printf("usage: <program name> <genome file> <generic options (see --help)\n");
    exit(1);
  }
  else
  {
    strcpy(genomefile,argv[1]);
   
    argv2=duplicateArgv(argc,argv,2);
    ReadPars(argc-1,argv2);// pass a selection of the command line to readpars
  
    char makedir[]="mkdir ";
    char makecommand[100];
    strcpy(makecommand, makedir);
    strcat(makecommand,despath);
  
    char testagent[100];
    
    if(stat(genomefile, &sb) != 0 || S_ISDIR(sb.st_mode) ){
      printf("genome file does not exist or is a directory. Exiting...\n");
      exit(1);
    }
    
    
    /// make directory if necessary ///
    if (stat(despath, &sb) == 0 && S_ISDIR(sb.st_mode)){//directory exists
      //do nothing;
    }
    
    else if(stat(despath, &sb) == 0 && !S_ISDIR(sb.st_mode)){//file exists but is not a directory
      printf("given destination is not a directory. give another name\n");
      exit(1);
    }
    else{ //directory does not exist: make it.
      if(system(makecommand)==-1){ //make directory for data
	printf("warning: could not make directory %s. Exiting...\n",despath);
	exit(1);
      }
    }
    
    ///make directory for movies if necessary ///
  /*  strcpy(extrapath,despath);
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
    }*/
    
  
  }
 
    
}

void doTransform(double **dat, int N, int Rep)
{
  //input: N real numbers, output: n/2+1
  fftw_complex *out; //output array
  fftw_plan p;
  double *in; //the array fed into the transform function
  
  /***array allocation ***/
  
  in=(double *)calloc(N, sizeof(double));

  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N/2+1));//output array
  
  /***********************/
  // the Plan: creates object containing all data to compute FFT
  p = fftw_plan_dft_r2c_1d(N,in, out, FFTW_MEASURE);
  //FFTW_FORWARD: from time to frequency (FFTW_BACKWARD does reverse)
  //FFTW_MEASURE: runs a number of FFT to find the best. FFTW_ESTIMATE just picks a reasonable one. For very big arrays, Estimate may be better but so far I haven't encountered serious delays
    
  //number of times we do the transform here (usually the number of genes)
  for(int i=0; i<Rep; i++){ 
    
    //transfer data from **dat to *in
    for(int j=0;j<N;j++)
    {
      in[j]=dat[i][j];
    }
    
    /*do the transform */
    fftw_execute(p); 
    
    store.push_back(vector < vector <double> >());
    
    /* store the output */
    for(int j=0;j<N/2+1;j++)
    {
      store[i].push_back(vector<double>());
      store[i][j].push_back(out[j][0]); //real part
      store[i][j].push_back(out[j][1]); //imaginary part
    }
  }
  
  
  fftw_destroy_plan(p);
  fftw_free(out);


}



int main(int argc, char **argv)
{
  FILE *f1, *f2, *f3, *f4;
  double sampling=1/(HT); //sampling rate
  char fname1[1000];
  char fname2[1000];
  
  Agent *A1;
  vector<Cell> C1, C2;
  vector< vector< vector<double> > > proteins, proteins2; //cell,  time, genes: stores time data (a cell doesn't do that)
  double **data;
  double **data2;
  
  /**read command line parameters */
  Start(argc,argv);
    
  /** initialize random numbers */
  fillgauss();
  dsfmt_init_gen_rand(&dsfmt,seedinitpop);//dsfmt_init_gen_rand(&dsfmt,11);
  
  /** allocate array for Fourier analysis **/
  data=(double **)calloc(NrGeneTypes, sizeof(double *));
  data[0]=(double *)calloc((NrDevSteps/*+EXTRATIME*/)*NrGeneTypes, sizeof(double));
  
  for(int i=1; i<NrGeneTypes; i++)
    data[i]=data[i-1]+(NrDevSteps/*+EXTRATIME*/);
  
  /************************************/
  /** start building the Agent, print basic info for good measure */
  A1=new Agent();
  A1->CreateAgentFromFile(genomefile,0); //easier than building the genome and network separately, have access to agent functions too
  A1->WriteSegmentEmbryology(0); //make more pictures
  /*A1->WriteBasicProperties(A1->agentid);
  A1->DetermineGenomeAndNetworkProperties();
  A1->WriteGenomeAndNetworkProperties(A1->agentid);
  A1->DetermineAttractorProperties();
  A1->WriteAttractorProperties(A1->agentid);
  A1->DetermineLoopAndMotifProperties();
  A1->WriteLoopAndMotifProperties(A1->agentid);
  
  A1->WriteNetworkProgression(0,extrapath,0);
  A1->WriteNetworkProgression(0,extrapath,10);   these functions just print the same as basic analysis, just clutter the directory*/
    
  
  /** Cells set 1: constant morphogen concentration, list with different concentrations **/
  /* initialize the cell list with artificial morphogens */
 
  for(int i=0; i<FIXCELLS; i++)
  {
    Cell c(i);
    c.SetCellState();
    c.proteinstates[0]=(double)(FIXCELLS-i)*100./(double)FIXCELLS;
    C1.push_back(c);
  }
  
  /* run the network */
   //open new file for data printing
   sprintf(fname2,"%s/DataFixedCells.dat",despath);
   f2=fopen(fname2,"w");
   
  for (int j=0;j<FIXCELLS;j++) //per cell, we run the network for the full duration
  {
    proteins.push_back(vector < vector <double > >()); //add a vector for this cell
    
        
    //set new cell to last state of previous cell
    if(j>0){
      for(int m=1;m<NrGeneTypes;m++) 
     	{
	  C1[j].proteinstates[m]=proteins[j-1][NrDevSteps-1][m];
      	}
    }
     /* run network for the cell */
    for(int i=0; i<NrDevSteps+EXTRATIME; i++)
    {
      proteins[j].push_back(vector <double>()); //add a vector of genes for this time point
      C1[j].UpdateCellState(A1->N); 
      fprintf(f2,"%d\t", i);
      
      for(int k=0; k<NrGeneTypes; k++){
	proteins[j][i].push_back(C1[j].proteinstates[k]); // store individual proteins at that time
	fprintf(f2,"%lf\t",C1[j].proteinstates[k]);
      }
      fprintf(f2,"\n");
    }
    
    fprintf(f2,"\n\n");
  }
  fclose(f2);
  /** ***************************** **/
  /** Cells set 2: decaying morphogen, different starting points **/
  /* initialize cell list with different starting points and decaying morphogens */
  sprintf(fname2,"%s/DataDecayingCells.dat",despath);
  f2=fopen(fname2,"w");
  for(int i=0; i<DECCELLS; i++)
  {
    Cell c(i);
    c.SetCellState();
    
    //set cell to one of the states of the cell with max morphogen concentration
    for(int m=0;m<NrGeneTypes;m++) 
    {
      c.proteinstates[m]=proteins[0][100+i*4][m];
    }
    C2.push_back(c);
  }
  
   /* run the network for the decaying cells */
  for (int j=0;j<DECCELLS;j++) //per cell, we run the network for the full duration
  {
    proteins2.push_back(vector < vector <double > >());
    
    for(int i=0; i<NrDevSteps; i++)
    {
      proteins2[j].push_back(vector <double>()); //add a vector of genes for this time point
      
      C2[j].UpdateCellState(A1->N);
      fprintf(f2,"%d\t", i);
      //update morphogen concentration
      C2[j].proteinstates[0]-=HT*morphdecay*C2[j].proteinstates[0];

      for(int k=0; k<NrGeneTypes; k++){
	fprintf(f2,"%lf\t",C2[j].proteinstates[k]);
	proteins2[j][i].push_back(C2[j].proteinstates[k]); //add individual proteins at that time
      }
      fprintf(f2,"\n");
    }
   fprintf(f2,"\n\n");
  }
  fclose(f2);
   
  /** print those cells that end up with high seg gene conc, and check when they reach equilibrium **/
  sprintf(fname2,"%s/HighCells.dat",despath);
  f2=fopen(fname2,"w");
  deque<int> success;
  int SSstep=0;//the timestep at which the simulations reach steady state
  int count=0;
  //printf("High cells: ");
  for (int i=0;i<50; i++)
  {
    if(proteins2[i][NrDevSteps-1][SegmGeneNr]-ThOn>-0.00001){
      success.push_back(i);
      fprintf(f2,"%d\n",i);//print the cells with high concentration
      while((proteins2[i][NrDevSteps-1][SegmGeneNr]-proteins2[i][NrDevSteps-1-count][SegmGeneNr])*
	(proteins2[i][NrDevSteps-1][SegmGeneNr]-proteins2[i][NrDevSteps-1-count][SegmGeneNr])<2. && count <NrDevSteps-1)
	count++;
      
      if(NrDevSteps-count>2)
	SSstep+=NrDevSteps-count;
      else
	success.pop_back();
      
      count=0;
    }
  }
  fclose(f2);
 //printf("SSstep: %d\n", SSstep);
  if(success.size())
    SSstep=(int)((double)SSstep/(double)success.size());//average time until reaching equilibrium
  else
    SSstep=NrDevSteps-1;
  
  sprintf(fname2,"%s/AverageFreezePoint.dat",despath);
  f2=fopen(fname2,"w"); 
  fprintf(f2,"%lf\n", proteins2[0][SSstep][0]);
  fclose(f2);
  //storage for the data of the decaying cells
  data2=(double **)calloc(NrGeneTypes, sizeof(double *));
  data2[0]=(double *)calloc(SSstep*NrGeneTypes, sizeof(double));
  
  for(int i=1; i<NrGeneTypes; i++)
    data2[i]=data2[i-1]+SSstep;
  
  /************************************/
  
  /** ***************************** **/
  
  /** ***************** **/
  /**  Fourier analysis **/
  /** ***************** **/
  
  /** Cells set 1 **/
  sprintf(fname1,"%s/FourierFixed.dat",despath);
  f1=fopen(fname1,"w");
  sprintf(fname1,"%s/FourierHeatmap.dat",despath);
  f2=fopen(fname1,"w");
  
  vector< vector< vector <double> > > corr; //cells series genes
  double prod;
  
  for (int i=0; i<FIXCELLS; i++){

      
    //append vector to corr for this cell
    corr.push_back(vector< vector< double> > ());
    
    //copy data into a temporary 2D array: for this "entire timeseries" analysis, we only use the series until NrDevSteps
    for(int j=0;j<NrGeneTypes; j++){
      for(int k=EXTRATIME; k<(NrDevSteps+EXTRATIME);k++){
	data[j][k-EXTRATIME]=proteins[i][k][j];
      }
    }
    
    doTransform(data, (NrDevSteps/*+EXTRATIME*/), NrGeneTypes); //fourier transform output is stored in the triple vector "store"
    
     //the series is printed in a column per gene
    for(int k=0; k<((NrDevSteps/*+EXTRATIME*/)/2+1);k++)
    {
      corr[i].push_back(vector< double> ());
      fprintf(f2,"%lf\t",C1[i].proteinstates[0]);
      fprintf(f1,"%lf\t",k*sampling/(NrDevSteps/*+EXTRATIME*/)); 
      fprintf(f2,"%lf\t",k*sampling/(NrDevSteps/*+EXTRATIME*/));
      //go through the genes
      for(int l=0; l<NrGeneTypes; l++) 
      {
	prod=sqrt(store[l][k][0]*store[l][k][0]+store[l][k][1]*store[l][k][1]);
	fprintf(f1,"%lf\t",prod); //print magnitude of the complex output
	fprintf(f2,"%lf\t", prod); //print magnitude of the complex output
	corr[i][k].push_back(prod); //store it for cross-correlations and such
      }
      fprintf(f1,"\n");
      fprintf(f2,"\n");
    }
    /*** do other stuff with the stored data? ***/
    
    store.clear();
    //store.resize();
    fprintf(f2,"\n");
    fprintf(f1,"\n\n");
}
  
 fclose(f2); 
 fclose(f1);
 
  
  
  
  /** Fourier joint time frequency (sliding window) analysis **/ //window of 300 steps, halfway overlap between slides: 11 points
  //double lastwindow[FIXCELLS][NrGeneTypes]; //stores for every cell and gene the amplitude of the last window. used to distinguish sustained and damped osc. for class.
  //construct a filename
  sprintf(fname1,"%s/FourierSlideWindow.dat",despath);
  f1=fopen(fname1,"w");
  for (int i=0; i<FIXCELLS; i++){
    
    for(int l=0; l<(NrDevSteps+EXTRATIME)/(WINDOWSIZE/2)-1; l++) //go through all window positions
    {
         
      //copy data into a temporary 2D array: for this "sliding window" analysis, we only use the series within the window
      for(int j=0;j<NrGeneTypes; j++){
	for(int k=0; k<WINDOWSIZE;k++){
	  data[j][k]=proteins[i][k+l*(WINDOWSIZE/2)][j];
	}
      }
      //fourier transform output is stored in the triple vector "store"
      doTransform(data, WINDOWSIZE, NrGeneTypes); 
      
      //the series is printed in a column per gene
    for(int k=1; k<(WINDOWSIZE/2+1);k++)
    {
           
      fprintf(f1,"%d\t",l); //time
      fprintf(f1,"%lf\t",k*sampling/(WINDOWSIZE)/* *2*3.14*/); //freq
      //go through the genes
      for(int l=0; l<NrGeneTypes; l++) 
      {
	fprintf(f1,"%lf\t", sqrt(store[l][k][0]*store[l][k][0]+store[l][k][1]*store[l][k][1]) ); //print magnitude of the complex output
	//if(k==(WINDOWSIZE/2))//we're at the last window: store the amplitude
	//{
	  //lastwindow[i][l]=sqrt(store[l][k][0]*store[l][k][0]+store[l][k][1]*store[l][k][1]);
	  //printf("lastwindow: %lf\n",lastwindow[i][l]);
	//}
      }
      fprintf(f1,"\n");
    }
    fprintf(f1,"\n");
    store.clear(); //empty storage for next window
    
    }
    fprintf(f1,"\n\n"); //give data for each cell in block
    
  }
  fclose(f1);
   
  
  /* find dominant frequencies for each gene and morphogen concentration */
  double maxim=0.;
  int freq;
  double topmax=0;
  int maxgene=0;
  
  sprintf(fname1,"%s/Morph_DomFreq.dat",despath);
  f1=fopen(fname1,"w");
  sprintf(fname2,"%s/Morph_Amp.dat",despath);
  f2=fopen(fname2,"w");
  sprintf(fname2,"%s/MaxOscGene.dat",despath); //stores for each gene its frequency and max amplitude in the growth zone
  f3=fopen(fname2,"w");
  sprintf(fname2,"%s/NativeGeneFreq.dat",despath); //stores for each gene its frequency and max amplitude in the growth zone
  f4=fopen(fname2,"w");
  
  for (int i=0; i<FIXCELLS; i++)
  {
    fprintf(f1, "%lf\t", C1[i].proteinstates[0]);
    fprintf(f2, "%lf\t", C1[i].proteinstates[0]);
        
    //find the peak for each gene
    for (int k=0; k<NrGeneTypes; k++)
    {
      for(int j=1; j<(NrDevSteps/*+EXTRATIME*/)/2+1; j++)//We don't start at the beginning: there is a ridiculous peak there
      {
	if((corr[i][j][k]-maxim)>1.){
	  maxim=corr[i][j][k];
	  freq=j;
	}
	
      }//timesteps loop
      
      if(i==0)//growth zone 
      {
	if(maxim>topmax) //find which gene oscillates the most
	{
	  topmax=maxim;
	  maxgene=k;
	}
	fprintf(f4, "%d %lf %lf\n",k,freq*sampling/(NrDevSteps/*+EXTRATIME*/), maxim);//print gene, frequency with max amp, max amplitude
      }
      if(maxim>1000) //we have sustained oscillations
      {
	fprintf(f1, "%lf\t", (double)freq*sampling/(NrDevSteps/*+EXTRATIME*/) /**2*3.14*/);
	fprintf(f2, "%lf\t", maxim/10000.);
      }
      else
      {
	fprintf(f1, "%lf\t", 0.);
	fprintf(f2, "%lf\t", 0.);
      }
      
      maxim=0.;
      freq=0;
    }//gene loop
    if(i==0)
    {
     fprintf(f3, "%d\n",maxgene); 
    }
    
    fprintf(f1,"\n");
    fprintf(f2,"\n");
   
  }//cell loop of different morphogen concentrations
  
  fclose(f1);
  fclose(f2);
  fclose(f3);
  fclose(f4);
  /** *********** **/
  corr.clear();
   
  // vector< vector< vector <double> > > corr; //cells series genes
  /** Cells set 2 **/
   //construct a filename
   sprintf(fname1,"%s/FourierDecay.dat",despath);
   f1=fopen(fname1,"w");
  for (int i=0; i<50; i++){
   
    
    //append vector to corr for this cell
    corr.push_back(vector< vector< double> > ());
    
    //copy data into a temporary 2D array
    for(int j=0;j<NrGeneTypes; j++){
      for(int k=0; k<SSstep;k++){
	data2[j][k]=proteins2[i][k][j];
      }
    }
    
    doTransform(data2, SSstep, NrGeneTypes); //fourier transform output is stored in the triple vector "store"
    
     //the series is printed in a column per gene
    for(int k=0; k<(SSstep/2+1);k++)
    {
      
      corr[i].push_back(vector< double> ());
      
      fprintf(f1,"%lf\t",k*sampling/SSstep /**2*3.14*/);
      //go through the genes
      for(int l=0; l<NrGeneTypes; l++) 
      {
	fprintf(f1,"%lf\t", sqrt(store[l][k][0]*store[l][k][0]+store[l][k][1]*store[l][k][1]) ); //print magnitude of the complex output
	corr[i][k].push_back(sqrt(store[l][k][0]*store[l][k][0]+store[l][k][1]*store[l][k][1])); //store it for cross-correlations and such
	
      }
      fprintf(f1,"\n");
    }
    /*** do other stuff with the stored data? ***/
    //check the profiles of genes between different cells
    store.clear();
    //store.resize();
    fprintf(f1,"\n\n");
   
  }
   fclose(f1);
  /** ************************ **/
  
   /** correlations and stuff **/
  
  vector < vector <double> > high (NrGeneTypes, vector <double> (SSstep/2+1, 0));
  vector < vector <double> > low (NrGeneTypes, vector <double> (SSstep/2+1, 0)); //genes, series
  int num=success.size();
  
  /*make average for high and low cells */
  for (int i=0; i<50; i++)
  {
    if(i==success[0]) //it is a high cell 
    {
      for (int p=0; p<NrGeneTypes; p++)
      {
	for(int w=0; w<SSstep/2+1;w++)
	{
	  high[p][w]+=corr[i][w][p]/(double)num;
	}
      }
      success.pop_front();
    }
    else
    {
      for (int p=0; p<NrGeneTypes; p++)
      {
	for(int w=0; w<SSstep/2+1;w++)
	{
	  low[p][w]+=corr[i][w][p]/(double)(50.-num);
	}
      }
    }
  }
  
  /**check for the max oscillating gene, how the frequency changes in the decaying cell **/
  //construct a filename
   sprintf(fname1,"%s/DecayPeakProfiles.dat",despath);
   f1=fopen(fname1,"w");
  
   //double locmax=0;
   double xposmax=0.;
   double prev=0;
   int found=0;
   
  for(int i=0; i<DECCELLS; i++)
  {
    fprintf(f1,"%d",i);
    for(int j=0;j<SSstep; j++)
    {
      if(proteins2[i][j][maxgene]>prev)
      {
	found=0;
	//locmax=proteins2[i][j][maxgene];
	xposmax=j;
      }
      
      else if(proteins2[i][j][maxgene]<prev && !found)
      {
	fprintf(f1,"\t%.2lf",xposmax);
	//locmax=0.0;
	found=1;
      }
      prev=proteins2[i][j][maxgene];
    }
    fprintf(f1,"\n");
  }
  fclose(f1);
  
  /*check for which genes, the difference is largest */
  
 /* sprintf(fname1,"%s/AverageHigh.dat",despath);
  f1=fopen(fname1,"w");
  
  sprintf(fname2,"%s/AverageLow.dat",despath);
  f2=fopen(fname2,"w");
  
  vector < vector <double> > diff (NrGeneTypes, vector <double> (SSstep/2+1, 0));
  double max=0.;
  int strong_gene=0;
  vector < double > overalldiff(NrGeneTypes,0);
  
  for(int w=0; w<SSstep/2+1;w++)
  {
    fprintf(f1, "%lf\t", w*sampling/SSstep // *2*3.14);
    fprintf(f2, "%lf\t", w*sampling/SSstep // *2*3.14);
    
    for (int p=0; p<NrGeneTypes; p++)
    {
      fprintf(f1, "%lf\t",high[p][w]);
      fprintf(f2, "%lf\t",low[p][w]);
      diff[p][w]=sqrt((high[p][w]-low[p][w])*(high[p][w]-low[p][w]));
      overalldiff[p]+=diff[p][w];
      
      if(diff[p][w]-max>0.001)
      {
	max=diff[p][w];
	strong_gene=p;
      }
    }
    fprintf(f1, "\n");
    fprintf(f2, "\n");
  }
  printf("gene with strongest difference peak: %d\n",strong_gene);
  for (int p=0; p<NrGeneTypes; p++)	
  {
    printf("gene %d: difference %lf\n", p, overalldiff[p]);
  }
  
  fclose(f1);
  fclose(f2);*/
  
 
 
  return 0;
  
}
