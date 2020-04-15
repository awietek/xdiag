#ifndef HYDRA_APPLICATIONS_HUBBARDFULLED_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDFULLED_OPTIONS_H_
#include <string>
#include "clara.hpp"


void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, int& nup, int& ndown, int& verbosity, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(couplingfile, "couplingfile")["-t"]["--couplingfile"]("name of couplingfile") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
    clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
    clara::Opt(verbosity, "verbosity")["-v"]["--verbosity"]("verbosity level, one of 0, 1 ,2 (optional, default 1)") |
    clara::Help(showhelp);


  auto cmd_args = parser.parse(clara::Args(argc,argv));
  if( !cmd_args ) 
    {
      std::cerr << "Error in command line: " << cmd_args.errorMessage() << std::endl;
      exit(EXIT_FAILURE);
    }
  else if (showhelp) 
    {
      parser.writeToStream(std::cout);
      exit(EXIT_SUCCESS);
    }
  else
    {
      std::cout <<
	"outfile     : " << outfile << std::endl <<
	"latticefile : " << latticefile << std::endl <<
	"couplingfile: " << couplingfile << std::endl <<
	"nup         : " << nup << std::endl <<
	"ndown       : " << ndown << std::endl <<
	"verbosity   : " << verbosity << std::endl <<
	"-----" << std::endl;
    }
}
#endif
