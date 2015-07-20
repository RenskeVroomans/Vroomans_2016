#ifndef GeneralHeader
#define GeneralHeader

using namespace std;
//using namespace __gnu_cxx;

#include <stdio.h>
#include <stdlib.h>
#include <ext/numeric>
#include <cmath>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <iterator>
#include <algorithm>
#include <boost/utility.hpp>
#include <cstring>
#include "dSFMT.h"

#define TRUE 1
#define FALSE 0

#define SUMINTEGRATION
//#define MULTIPLYINTEGRATION

#define RUN
//#define RERUN
//#define NOANCESTRY
//#define INITGENOME

//#define INITGRAD
#define TWOMORPHS

#define GRADIENT
//#define FREEMORPH
//#define POSTERIORSIGNAL
//#define MORPHDIFF
//#define WAVEFRONT


extern char despath[500];

/** vars that have to be constant **/
//in Population.cc
const int NI=30;
const int NJ=30;
//used in Agent.cc
const int NrGeneTypes=16;
#ifdef TWOMORPHS
const int NrMatGeneTypes=2;//gene 0
#else
const int NrMatGeneTypes=1;//gene 0
#endif

const int NrSignGeneTypes=0;//4;//2;// genes that diffuse or signal to neighbours
const int GrowGeneNr=4;//NrMatGeneTypes+NrSignGeneTypes; //the gene that regulates growth
const int SegmGeneNr=5;//NrMatGeneTypes+NrSignGeneTypes;//gene 5 (was 7) the "output" gene 

const int NrStorages=121;//1+(NrDevSteps/StorageInt);//not read from parfile
const int NrFinalCells=150;//134;//125;//200
const int minbandsize=7;//9;//8;
const int MaxNrBands=NrFinalCells/minbandsize; //not read from parfile

/** can be read from parfile **/
//in Population.cc
extern int InitNrAg; //50;
extern double deathrate;

//in Agent.cc
extern int StorageInt;
extern int NrDevSteps;

extern double ThOn; //threshold value for a gene to be considered 'on' (in fitness/stripe calculation
extern double ThOff;//threshold val for gene to be considered 'off'
extern double Emax;
extern double Decay;
extern double morphdecay; //0.05
extern double DifCoef;
extern int Nrdiffsteps;
extern int divinterval;//double divprob=0.83;
extern double HT;//0.2;//1.;//0.5;
extern int MaintInt;//10;//20;
extern int InitNrCells;//128;;//1;
extern int targetsize;
extern int storecellspec; //whether you want cellid specific time info
extern int averagepattern; //whether the fitness should be measured over an average of the last MaintInt timesteps
extern double genepen;
extern double tfbspen;
extern double sizebonus;
extern double sizepen;
extern double stablesizepen;
extern double regpen;

//used in Genome.cc
extern int AvCon;//5;//3; //average connectivity of a node (?)
extern double fac;
extern double fac2;

//used in Network.cc
extern double H;//30;//60.;
extern int N;//4;//2;


/***not read from parfile, but modified upon reading***/
extern double probnontandgenedupl;
extern double probgenedel;//0.00120*fac;//0.000090*fac;
//double probidimpswitch=0.000020*fac;
extern double probtfbsweightrev;
extern double probtfbstypeswitch;
extern double probnontandtfbsdupl;//0.000020*fac;//per tfbs
extern double probtfbsinnov;//per genome
extern double probtfbsdel;//0.0000375*fac;//0.000040*fac
/**************************/

//in World.cc
extern int seedinitpop;
extern int seedmutations;
extern int seedoutput;
extern int NrGenerations;//10000;

//rest
extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }

extern dsfmt_t dsfmt_output;
inline double uniform_output() { return dsfmt_genrand_close_open(&dsfmt_output); }

extern double gauss[10000];

#endif
