#include "Header.hh"
#include "Population.hh"
#include "Agent.hh"
#include "Genome.hh"
#include "Gene.hh"
#include "TFBS.hh"
#include "Misc.hh"
#include "dSFMT.h"
#include "IO.hh"

#include <sys/stat.h>

dsfmt_t dsfmt;
dsfmt_t dsfmt_output;


double gauss[10000];



void Start(int argc,char **argv)
{
  struct stat sb;        //allows to check for the existence of directories
  
  /*if(argc<3)
    {
      printf("usage: world destinationpath seed\n");
      exit(1);
    }
  else
    {
      strcpy(despath,argv[1]);
      seedinitpop=atoi(argv[2]);
    }*/
  ReadPars(argc,argv);
  
  char makedir[]="mkdir ";
  char makecommand[50];
  strcpy(makecommand, makedir);
  strcat(makecommand,despath);
  
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

}



int main(int argc, char **argv)
{
  int i;
  Population *P;
 // char ancestryfile[800]="PopAncestry10000";
  //char genomeinputfile[800]="CodedGenomeAgent0000268555";
  char command[800];
  char direc[800];

  fillgauss();
  Start(argc,argv);

#ifdef RUN
  P=new Population();
  dsfmt_init_gen_rand(&dsfmt, seedinitpop);
  (*P).InitPopulation();
  dsfmt_init_gen_rand(&dsfmt, seedmutations);
  printf("all is initialised\n");
  (*P).WritePopFitnessToFile(0);
  (*P).WriteGenomeLengthToFile(0);
  (*P).WriteBandsToFile(0);
  (*P).WriteFittestMatricesToFile(0);
#endif
#ifdef INITGENOME
  P=new Population();
  dsfmt_init_gen_rand(&dsfmt, seedinitpop);
  (*P).InitPopulationFromGenome(genomeinputfile);
  dsfmt_init_gen_rand(&dsfmt, seedmutations);
  printf("all is initialised\n");
  (*P).WritePopFitnessToFile(0);
  (*P).WriteGenomeLengthToFile(0);
  (*P).WriteBandsToFile(0);
  (*P).WriteFittestMatricesToFile(0);
#endif
#ifdef RERUN
  P=new Population();
  sprintf(direc,"AncestryGenomes");
  sprintf(command, "mkdir %s/%s",despath,direc);
  if(system(command)==-1)
  {
    printf("AncestryGenomes: warning: could not make dir. Exiting...\n");
    exit(1);
  }
  (*P).ReadAncestorsFromFile(ancestryfile);
  dsfmt_init_gen_rand(&dsfmt, seedinitpop);
  (*P).InitPopulation();
  dsfmt_init_gen_rand(&dsfmt, seedmutations);
  printf("all is initialised\n");
  (*P).WriteBandsToFile(0);
  (*P).WriteFittestMatricesToFile(0);
#endif
#ifdef NOANCESTRY
  P=new Population();
  dsfmt_init_gen_rand(&dsfmt, seedinitpop);
  (*P).InitPopulation();
  dsfmt_init_gen_rand(&dsfmt, seedmutations);
  printf("all is initialised\n");
  (*P).WritePopFitnessToFile(0);
  (*P).WriteGenomeLengthToFile(0);
  (*P).WriteBandsToFile(0);
  (*P).WriteFittestMatricesToFile(0);
#endif

  dsfmt_init_gen_rand(&dsfmt_output,seedoutput);
  
  for(i=1;i<NrGenerations;i++)
    {
      (*P).Step(i);  
   
#ifdef RUN
      if(i%25==0)
	{
	  (*P).WritePopFitnessToFile(i);
	  (*P).WriteGenomeLengthToFile(i);
	  (*P).WriteBandsToFile(i);
	}
      if(i%100==0)
	{
	  (*P).WriteFittestMatricesToFile(i);
	}
      if(i%1000==0)
      {
	(*P).WriteAncestriesToFile(i);
	(*P).WriteFieldToFile(i);
      }
#endif
#ifdef RERUN
      if(i%10==0)
	(*P).WriteBandsToFile(i);
      if(i%100==0)
	(*P).WriteFittestMatricesToFile(i);
#endif
#ifdef NOANCESTRY
      if(i%10==0)
	 {
	  (*P).WritePopFitnessToFile(i);
	  (*P).WriteGenomeLengthToFile(i);
	  (*P).WriteBandsToFile(i);
	}
      if(i%100==0)
	{
	  (*P).WriteFittestMatricesToFile(i);
	}
#endif
    }

  //delete P;

  return(0);
}
