#include <stdio.h>
#include <stdlib.h>
#include "Header.hh"

/** can be read from parfile **/
char despath[500];
#ifdef RERUN
char ancestryfile[500];
#endif
//in Population.cc
int InitNrAg=1; //50;
double deathrate=0.5;

//in Agent.cc
double ThOn=80.; //threshold value for a gene to be considered 'on' (in fitness/stripe calculation
double ThOff=20.;//threshold val for gene to be considered 'off'
double morphdecay=0.025; //0.05
double DifCoef=0.;
int Nrdiffsteps=10;
int divinterval=5;//double divprob=0.83;
double HT=0.2;//0.2;//1.;//0.5;
int MaintInt=20;//10;//20;
int InitNrCells=5;//128;;//1;
int storecellspec=1; //whether you want cellid specific time info
int averagepattern=0;
int NrDevSteps=600;//600;//120;//150;
int targetsize=150;
int StorageInt=5;//5;//1;//6;

double genepen=0.00001; //penalty for nr of genes
double tfbspen=0.000001;//penalty for nr of tfbs
double sizebonus=0.05;//bonus for growing towards target size
double sizepen=2; //penalty for exceeding target size
double stablesizepen=0.;
double regpen=0.;

//used in Genome.cc
int AvCon=3;//5;//3; //average connectivity of a node (?)
double fac=100;
double fac2=1.;
double Dmin=0.05;
double Dmax=0.9;
double Hmin=10;
double Hmax=400;

//used in Network.cc
double H=60.;//30;//60.;
int N=2;//4;//2;
double Estart=70.0;
double Decay=0.3;
double Emin=10.;
double Emax=120.;

/***not read from parfile, but modified upon reading***/
double probnontandgenedupl=0.000060*fac;
double probgenedel=0.000090*fac;//0.00120*fac;//0.000090*fac;
double probgeneDDchange=0.000030*fac;
double probgeneEEchange=0.000030*fac;
//double probidimpswitch=0.000020*fac;
double probtfbsweightrev=0.000010*fac*fac2;
double probtfbstypeswitch=0.000010*fac*fac2;
double probnontandtfbsdupl=0.000015*fac*fac2;//0.000020*fac;//per tfbs
double probtfbsinnov=0.000010*fac*fac2;//per genome
double probtfbsdel=0.000040*fac*fac2;//0.0000375*fac;//0.000040*fac
double probtfbsHHchange=0.000020*fac*fac2;
/**************************/

//in World.cc
int seedinitpop=3;
int seedmutations=9;
int seedoutput=11;
int NrGenerations=1;//10000;
