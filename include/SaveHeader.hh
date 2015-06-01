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


#define RUN
//#define RERUN
//#define NOANCESTRY
//#define INITGENOME


#define GRADIENT
//#define WAVEFRONT

extern char despath[500];

/** vars that have to be constant **/
//in Population.cc
const int NI=5;//30;
const int NJ=5;//30;
//used in Agent.cc
const int NrGeneTypes=16;
const int NrMatGeneTypes=1;//gene 0
const int NrSignGeneTypes=0;//4;//2;// genes that diffuse or signal to neighbours
const int SegmGeneNr=5;//NrMatGeneTypes+NrSignGeneTypes;//gene 5 (was 7) the "output" gene 
const int NrDevSteps=600;//600;//120;//150;
const int StorageInt=5;//5;//1;//6;
const int NrStorages=1+(NrDevSteps/StorageInt);//not read from parfile
const int NrFinalCells=125;//128;
const int minbandsize=7;//9;//8;
const int MaxNrBands=NrFinalCells/minbandsize; //not read from parfile

/** can be read from parfile **/
//in Population.cc
const int InitNrAg=1; //50;
const double deathrate=0.5;

//in Agent.cc
const double ThOn=80.; //threshold value for a gene to be considered 'on' (in fitness/stripe calculation
const double ThOff=20.;//threshold val for gene to be considered 'off'
const double Emax=100.0;
const double Decay=0.3;
//double morphprod=0.5;
const double morphdecay=0.025; //0.05
const int growzonesize=5;
const int divinterval=5;//double divprob=0.83;
const double HT=0.2;//0.2;//1.;//0.5;
const int MaintInt=20;//10;//20;
const int InitNrCells=5;//128;;//1;
const int storecellspec=1; //whether you want cellid specific time info

//used in Genome.cc
const int AvCon=3;//5;//3; //average connectivity of a node (?)
const double fac=100;
const double fac2=1.;

//used in Network.cc
const double H=60.;//30;//60.;
const int N=2;//4;//2;


/***not read from parfile, but modified upon reading***/
const double probnontandgenedupl=0.000060*fac;
const double probgenedel=0.000090*fac;//0.00120*fac;//0.000090*fac;
//double probidimpswitch=0.000020*fac;
const double probtfbsweightrev=0.000010*fac*fac2;
const double probtfbstypeswitch=0.000010*fac*fac2;
const double probnontandtfbsdupl=0.000015*fac*fac2;//0.000020*fac;//per tfbs
const double probtfbsinnov=0.000010*fac*fac2;//per genome
const double probtfbsdel=0.000040*fac*fac2;//0.0000375*fac;//0.000040*fac
/**************************/

//in World.cc

const int seedmutations=9;
const int seedoutput=11;
const int NrGenerations=1;//10000;

//rest
extern dsfmt_t dsfmt;
inline double uniform() { return dsfmt_genrand_close_open(&dsfmt); }

extern dsfmt_t dsfmt_output;
inline double uniform_output() { return dsfmt_genrand_close_open(&dsfmt_output); }

extern double gauss[10000];

#endif
