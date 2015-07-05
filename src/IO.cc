//#ifndef IOHEADER_GUARD  // Top of file
//#define IOHEADER_GUARD

#include "IO.hh"
//#include <stdio.h>
//#include <stdlib.h>
#include <boost/program_options.hpp>
#include <iostream> 
#include <string>
#include <vector>
#include <fstream>

using namespace boost;
namespace po = boost::program_options;
using namespace std;

char *convert(const std::string & s)
{
   char *pc = new char[s.size()+1];
   std::strcpy(pc, s.c_str());
   return pc; 
}


int ReadPars(int argc, char* argv[])
{
  /** ******command line reading******** **/
  
  // Declare the supported options.
  po::options_description generic("command line options");
  //this function actually lists the options to read from command line
  generic.add_options()
  ("help", "produce help message")
  ("parfile",po::value< vector <string> >(),"file with the parameters")
  ("despath,d",po::value< vector <string> >()->required(),"destination path")
  ("seedinitpop,s", po::value<int>(&seedinitpop)->default_value(503))
  ;
  
  
  //to make sure that you don't need to say --parfile="name.cfg", but just name.cfg
  po::positional_options_description p;
  p.add("parfile", -1);
  
  po::variables_map vm;//stores values of options, and can store values of arbitrary types. 
  //store, parse_command_line and notify functions cause vm to contain all the options found on the command line.
  po::store(po::command_line_parser(argc, argv).options(generic).positional(p).run(), vm);
  
  if (vm.count("help")) { //checks if option was specified
    cout << generic << "\n";
    //return 1;
    exit(1);
  }
  vector <string> vs=vm["despath"].as< vector<string> >();
  vector<char*>  vc;
  
  transform(vs.begin(), vs.end(), std::back_inserter(vc), convert);
  strcpy(despath,&vs[0][0]);
  cout <<"put files in folder: "<<&vs[0][0]<<endl;
  
  if(!(vm.count("parfile"))){
    cout << "no parfile specified. Using default values at own risk..." <<endl;
    //exit(1);
  }
 
  /** ******parfile reading******** **/
  
  else{
    
    /**specify the parameters here!!**/
    po::options_description config("parameter values");
    config.add_options()
    //("NI", po::value<int>(&NI)->default_value(30), "nr of rows")
    //("NJ", po::value<int>(&NJ)->default_value(30), "nr of cols")
    ("InitNrAg", po::value<int>(&InitNrAg)->default_value(50), "nr of agents to start with")
    ("deathrate", po::value<double>(&deathrate)->default_value(0.5), "probability of agent dying")
    //used in Agent.cc
    //("NrGeneTypes", po::value<int>(&NrGeneTypes)->default_value(16), "nr gene types")
    //("NrMatGeneTypes", po::value<int>(&NrMatGeneTypes)->default_value(1), "nr maternal genes")
    //("NrSignGeneTypes", po::value<int>(&NrSignGeneTypes)->default_value(0), "nr gene types that diffuse or signal to neighbours")
    //("SegmGeneNr", po::value<int>(&SegmGeneNr)->default_value(5), "gene responsible for segmentation")//NrMatGeneTypes+NrSignGeneTypes
    ("ThOn", po::value<double>(&ThOn)->default_value(80.), "threshold val for gene to be considered on (in DetermineFitness)")
    ("ThOff", po::value<double>(&ThOff)->default_value(20.), "threshold val for gene to be considered off (in DetermineFitness)")
    ("Estart", po::value<double>(&Estart)->default_value(70.))
    ("Decay", po::value<double>(&Decay)->default_value(0.3))
    ("morphdecay", po::value<double>(&morphdecay)->default_value(0.05))
    ("DifCoef", po::value<double>(&DifCoef)->default_value(0.0))
    ("Nrdiffsteps", po::value<int>(&Nrdiffsteps)->default_value(10))
    ("InitNrCells", po::value<int>(&InitNrCells)->default_value(5), "nr of cells to start with")
    ("divinterval", po::value<int>(&divinterval)->default_value(5), "nr of steps between divisions of posterior cell")
    ("HT", po::value<double>(&HT)->default_value(0.2),"time step of integration")
    ("NrDevSteps", po::value<int>(&NrDevSteps)->default_value(600), "nr of steps to develop the agent")
    //("StorageInt", po::value<int>(&StorageInt)->default_value(5), "nr of steps between storage of the protein levels")
    ("MaintInt", po::value<int>(&MaintInt)->default_value(20),"distance from final timestep to store pattern for stability check")
    ("targetsize", po::value<int>(&targetsize)->default_value(100),"target nr of cells to which the animal should grow")
    //("minbandsize", po::value<int>(&minbandsize)->default_value(7),"min size of proper segment")
    ("storecellspec", po::value<int>(&storecellspec)->default_value(1),"whether time data should be printed for specific cell ids")
    ("averagepattern", po::value<int>(&averagepattern)->default_value(0),"whether the fitness should be measured over an average of the last MaintInt timesteps")
    ("genepen", po::value<double>(&genepen)->default_value(0.00001),"fitness penalty for nr of genes")
    ("tfbspen", po::value<double>(&tfbspen)->default_value(0.000001),"fitness penalty for nr of tfbs")
    ("sizebonus", po::value<double>(&sizebonus)->default_value(0.1),"fitness bonus for being closer to target size")
    ("sizepen", po::value<double>(&sizepen)->default_value(2.),"fitness penalty for exceeding target size")
    ("stablesizepen", po::value<double>(&stablesizepen)->default_value(0.),"fitness penalty for dividing within maintenance period")
    //used in Genome.cc
    ("AvCon", po::value<int>(&AvCon)->default_value(3), "average connectivity of a node?")
    ("fac", po::value<double>(&fac)->default_value(100))
    ("fac2", po::value<double>(&fac2)->default_value(1.))
    //used in Network.cc
    ("H", po::value<double>(&H)->default_value(60.))
    ("N", po::value<int>(&N)->default_value(2))
    ("Dmin", po::value<double>(&Dmin)->default_value(0.05))
    ("Dmax", po::value<double>(&Dmax)->default_value(0.9))
    ("Hmin", po::value<double>(&Dmin)->default_value(10.))
    ("Hmax", po::value<double>(&Dmax)->default_value(400.))
    ("Emin", po::value<double>(&Emin)->default_value(10.))
    ("Emax", po::value<double>(&Emax)->default_value(120.))
    //in World.cc
    ("seedmutations", po::value<int>(&seedmutations)->default_value(9))
    ("seedoutput", po::value<int>(&seedoutput)->default_value(11))
    ("NrGenerations", po::value<int>(&NrGenerations)->default_value(10000))
    
    
    ;

    /******required type conversions: from vector<string> to ifstream for parse_config_file()*******/
    vector <string> vs=vm["parfile"].as< vector<string> >();
    vector<char*>  vc;
    transform(vs.begin(), vs.end(), std::back_inserter(vc), convert);
  
    cout<<"parfile used: "<<&vc[0][0]<<endl;
    ifstream parf(&vc[0][0]);
    /*************/
  
    po::store(po::parse_config_file(parf, config), vm);
    
    po::notify(vm);    
    
    /***variables with calculations***/
    
    //in Agent.cc:
    StorageInt=(int)ceil((double)NrDevSteps/(double)(NrStorages-1));
 
    //cout << "storageint: " << StorageInt << ", NrDevSteps: " << NrDevSteps << endl;
 
    probnontandgenedupl=0.000060*fac;
    probgenedel=0.000090*fac;//0.00120*fac;//0.000090*fac;
    probgeneDDchange=0.000030*fac;
    probgeneEEchange=0.000030*fac;
    //const double probidimpswitch=0.000020*fac;
    probtfbsweightrev=0.000010*fac*fac2;
    probtfbstypeswitch=0.000010*fac*fac2;
    probnontandtfbsdupl=0.000015*fac*fac2;//0.000020*fac;//per tfbs
    probtfbsinnov=0.000010*fac*fac2;//per genome
    probtfbsdel=0.000040*fac*fac2;//0.0000375*fac;//0.000040*fac
    probtfbsHHchange=0.000020*fac*fac2;
  
    
    
  }
  
  return 0;

}
//#endif