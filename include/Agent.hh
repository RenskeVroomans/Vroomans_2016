#ifndef AgentHeader
#define AgentHeader

#include "Header.hh"
#include "Misc.hh"
#include "Genome.hh"
#include "Network.hh"
#include "Cell.hh"
#include "stdlib.h"

class Agent
{
 public:
  Genome *G;//1 per agent!
  Network *N;//1 per agent!
  list<Cell> cells;
  int E[NrStorages][NrFinalCells][NrGeneTypes]; //store spatiotemporal expression pattern
  double averseg[NrFinalCells];//to store the average expression of the segmentation gene per position
  double varseg[NrFinalCells];
  int cellpattern[NrStorages][NrFinalCells/10][NrGeneTypes]; //stores pattern per cell, not per position
  int types[NrStorages][NrFinalCells];//store spatiotemporal celltype pattern
  int ages[NrStorages][NrFinalCells];//store cell ages
  int bands[NrFinalCells][4];
  int bandtypes[NrFinalCells];
  int nrbands;
  int nrlongbands;
  int agentid;
  int parentid;
  int tbirth;
  double nonexpfitness;
  double fitness;
  double glpenalty;
  double instpenalty;
  double shortsegpenalty;
  double sizefit;
  double regpenalty;
  int anrcells_;
  int maintsize;
  double wfspeed;
  //double divprob; //individual division probability for stochastic divisions

  list<int> ancestry;
  typedef std::list<int>::iterator iter;
  
  Agent();
  ~Agent();
  void CreateAgentFromScratch(int agid,int t);
  void CreateAgentFromGenome(int agid,int t,char *fname);
  void CloneAgentFromAgent(Agent *A,int agid,int t);
  void CreateAgentFromParent(Agent *A,int agid,int t);
  void InitAgent(int aid,int pid,int t);
  void DevelopAgent();
  void FormZygote();
  void RemoveGradient();
  void MorphogenDiffusion(int prot);
  void CellCellSignalling(int t);
  void IntracellularDynamics();
  void DivideCells();
  void StoreAgentState(int I);
  void StoreCellStates(int I);
  void MaintenanceIntracellularDynamics(int i);
  void DetermineFitness(int mode);
  void WriteGenome(char *dirname);
  int CellTypeToColor(int type);
  int SegmentToColor(int type, int max);
  int ProteinStateToCellType(double proteinstates[NrGeneTypes]);
  void WriteEmbryology(char *dirname);
  void WriteDivisionProfile(char *dirname, int c);
  void WriteSignalProfiles(char *dirname);
  void WriteFitnessDetails(char *dirname);
  
  
  //for debugging
  void UpdateCellAges();
  void PrintAgentData(int tijd);
  void PrintCellAges();
  
  //for analysis: In AgentAnalysis.cc
  void CreateAgentFromFile(char *despath,int agentidnr);
  void InitAgentFromFile(char *despath,int agentidnr);
  void FinishAgentFromFile();
  void DetermineGenotypeSimilarity(Agent *M,int *genosim);
  void CloneAgentFromAgentForSafety(Agent *A,int agentidnr,int flag);
  //int DifferentExpressionPattern(double p1[NrGeneTypes],double p2[NrGeneTypes]);
  void FindProteinOscillations(int *array);
  int DifferentEquilibrium(double p1[NrGeneTypes],double p2[NrGeneTypes]);
  int CompareAgents(Agent *A, int c);
  void WriteGenome(int c);
  void WriteNetwork(char* c);
  void WriteNetworkProgression(int c, char* destdir, int cell);
  //void WriteMovieEmbryology();
  void WriteFullAgeProfile(char *subdir, int c);
  void WriteTimepointAgeProfile(char *subdir, int c, int timepoint);
  void WriteCloneProfile(char *dirname, int c, int timepoint);
  void WriteCloneTimecourse(char *subdir,int c, int timepoint);
  void WriteGeneEmbryology(char *subdir, int c, int gene);
  void WriteTimepointGene(char *subdir, int c, int gene, int timepoint);
  void WriteVideo(int gene);
  void WriteVideoPicture(int gene, int time, int max);
  void WriteExtraSignalProfiles(int c);
  void WriteBasicProperties(int aid);
  void DetermineGenomeAndNetworkProperties(char *name);
  void removeDoubles_type();
  void removeDoubles_id(int finalcall);
  void WriteGenomeAndNetworkProperties(int aid);
  void DetermineAttractorProperties();
  void WriteAttractorProperties(int aid);
  void DetermineLoopAndMotifProperties();
  void WriteLoopAndMotifProperties(int aid, char *name);
  void DetermineGenoNetworkPhenoDistance(Agent *M,double *genodist,double *phenodist);
  void WriteGenoNetworkPhenoDistance(int aid,double genodist,double phenodist,double cumgenodist,double cumphenodist);
  
  int nrlonganddifbands;//domain needs to express at least 1 id gene, but minimum is 1 domain even if nowhere id gene
  int nrhomogbands;
  int nrdifbutrepeatedbands;
  int nrequilibria;
  double equilibria[1000][NrGeneTypes];
  double sizeequilibria[1000];
  int nrequsedtotal,nrequsedend;
  double avgeqsize;
  int nrplustfbs,nrmintfbs;
  double avgindegree,maxindegree,minindegree,stdindegree;
  double avgoutdegree,maxoutdegree,minoutdegree,stdoutdegree;
  double avgtotdegree,maxtotdegree,mintotdegree,stdtotdegree;
  int nrposauto, nrnegauto, nrpos, nrneg;
  int compreg;
  list< vector<int> > loopstore, counterstore; //store loops (gene type or gene ID)
  vector <int> tallying; 
  
};
#endif

