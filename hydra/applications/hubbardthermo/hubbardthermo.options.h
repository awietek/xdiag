#ifndef HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_

#include <string>
#include "clara.hpp"

void parse_cmdline(double& T, double& t, double& U, std::string& outfile, 
		   std::string& latticefile, std::string& temperaturefile, 
		   int& nup, int& ndown, int& np, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(T, "temperature")["-t"]["--temperature"]("temperature (ignored if temperaturefile is given)") |
    clara::Opt(t, "t")["-t"]["--t"]("Hubbard t") |
    clara::Opt(U, "U")["-U"]["--U"]("Hubbard U") |
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(temperaturefile, "temperaturefile")["-t"]["--temperaturefile"]("name of temperaturefile") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
    clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
    clara::Opt(np, "np")["-p"]["--np"]("number of particles (optional, default half filling)") |
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
	"t              : " << t << std::endl <<
	"U              : " << U << std::endl <<
	"outfile        : " << outfile << std::endl <<
	"latticefile    : " << latticefile << std::endl <<
	"temperaturefile: " << temperaturefile << std::endl <<
	"nup            : " << nup << std::endl <<
	"ndown          : " << ndown << std::endl <<
	"np             : " << np << std::endl <<
	"-----" << std::endl;
    }
}

#endif
