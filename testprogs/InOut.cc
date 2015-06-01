//#include "header.hh"
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

int main(int argc, char* argv[])
{
  
  int var1, var2;
  
  // Declare the supported options.
  po::options_description generic("Allowed options");
  //this function actually lists the options to read
  generic.add_options()
  ("help", "produce help message")
  ("parfile",po::value< vector <string> >(),"file with the parameters")
  ;
  po::positional_options_description p;
  p.add("parfile", -1);
  
  po::variables_map vm;//stores values of options, and can store values of arbitrary types. 
  //store, parse_command_line and notify functions cause vm to contain all the options found on the command line.
  po::store(po::command_line_parser(argc, argv).options(generic).positional(p).run(), vm);
  
  if (vm.count("help")) { //checks if option was specified
    cout << generic << "\n";
    return 1;
  }
    
  if(!(vm.count("parfile"))){
    cout << "please specify a parfile" <<endl;
    exit(1);
  }
    
  po::options_description config("Configuration");
  config.add_options()
  ("randvar", po::value<int>(&var1)->default_value(10), "just a var")
  ("randvar2", po::value<int>(&var2)->required(), "just a var")
  ;
  
  
  
  
  /*************/
  vector <string> vs=vm["parfile"].as< vector<string> >();
  vector<char*>  vc;
  transform(vs.begin(), vs.end(), std::back_inserter(vc), convert);
  
  cout<<&vc[0][0]<<endl;
   ifstream parf(&vc[0][0]);
  //parf.open(parf,ifstream::in);//&vc[0][0]
  /*************/
  
  po::store(po::parse_config_file(parf, config), vm);
  
  po::notify(vm);    
  
  //cout << "Input files are: "<< vs << "\n"; 
	 
  
  cout << var1 << var2 << endl;
  return 0;

}
