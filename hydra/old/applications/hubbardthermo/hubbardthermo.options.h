#ifndef HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_

#include <string>
#include "clara.hpp"

void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, std::string& ensemble, std::string& temperaturefile, int& np, double& mu, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(couplingfile, "couplingfile")["-c"]["--couplingfile"]("name of couplingfile (default: latticefile)") |
    clara::Opt(ensemble, "ensemble")["-e"]["--ensemble"]("ensemble to use, either canonical or grandcanonical (default: canonical)") |
    clara::Opt(temperaturefile, "temperaturefile")["-t"]["--temperaturefile"]("name of temperaturefile") |
    clara::Opt(np, "np")["-p"]["--np"]("number of particles for canonical ensemble (optional, default half filling)") |
    clara::Opt(mu, "mu")["-m"]["--mu"]("chemical potential for grandcanonical ensemble (optional, default 0.)") |
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
	"outfile        : " << outfile << std::endl <<
	"latticefile    : " << latticefile << std::endl <<
	"couplingfile   : " << couplingfile << std::endl <<
	"ensemble       : " << ensemble << std::endl <<
	"temperaturefile: " << temperaturefile << std::endl <<
	"np             : " << np << std::endl <<
	"mu             : " << mu << std::endl <<
	"-----" << std::endl;
    }
}
#endif
