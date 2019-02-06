#ifndef HYDRA_APPLICATIONS_HEISENBERGTHERMO_OPTIONS_H_
#define HYDRA_APPLICATIONS_HEISENBERGTHERMO_OPTIONS_H_

#include <string>
#include "clara.hpp"

void parse_cmdline(double& T, double& J, std::string& outfile, 
		   std::string& latticefile, std::string& temperaturefile, 
		   int& nup, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(T, "temperature")["-t"]["--temperature"]("temperature (ignored if temperaturefile is given)") |
    clara::Opt(J, "t")["-t"]["--t"]("coupling constant J") |
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(temperaturefile, "temperaturefile")["-t"]["--temperaturefile"]("name of temperaturefile") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up spins (optional, default: all sectors are combined)") |
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
	"T              : " << T << std::endl <<
	"J              : " << J << std::endl <<
	"outfile        : " << outfile << std::endl <<
	"latticefile    : " << latticefile << std::endl <<
	"temperaturefile: " << temperaturefile << std::endl <<
	"nup            : " << nup << std::endl <<
	"-----" << std::endl;
    }
}
#endif
