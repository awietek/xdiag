#ifndef HYDRA_APPLICATIONS_HEISENBERGED_OPTIONS_H_
#define HYDRA_APPLICATIONS_HEISENBERGED_OPTIONS_H_

#include <string>
#include "clara.hpp"

void parse_cmdline(double& J, std::string& outfile, std::string& latticefile, 
		   int& nup, std::string& representation, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(J, "J")["-J"]["--J"]("hopping constant t") |
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of upspins (optional, default: half filled)") |
    clara::Opt(representation, "representation")["-r"]["--representation"]("space group representation (optional: default no space group symmetries applied)") |
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
	"J             : " << J << std::endl <<
	"outfile       : " << outfile << std::endl <<
	"latticefile   : " << latticefile << std::endl <<
	"nup           : " << nup << std::endl <<
	"representation: " << representation << std::endl <<
	"-----" << std::endl;
    }
}

#endif
