#ifndef HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_

#include <string>
#include <sstream>

#include <lila/allmpi.h>
#include "clara.hpp"

void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, 
		   std::string& ensemble, std::string& temperaturefile, int& np, double& mu, 
		   int& seed, double& mintemperature, double& precision, int& iters, 
		   int& verbosity, int& argc, char** argv)

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
    clara::Opt(seed, "seed")["-s"]["--seed"]("random seed (optional, default 1)") |
    clara::Opt(mintemperature, "mintemperature")["-m"]["--mintemperature"]("minumum temperature for convergence (optional, default 0.01)") |
    clara::Opt(precision, "precision")["-p"]["--precision"]("precision of Lanczos imag. time evolution (optional, default 1e-12)") |
    clara::Opt(iters, "iters")["-i"]["--iters"]("maximum number of Lanczos iterations performed (optional, default 1000)") |
    clara::Opt(verbosity, "verbosity")["-v"]["--verbosity"]("verbosity level, one of 0, 1 ,2 (optional, default 1)") |
    clara::Help(showhelp);



  auto cmd_args = parser.parse(clara::Args(argc,argv));

  lila::loggermpi.set_verbosity(verbosity);
  if( !cmd_args ) 
    {
      lila::loggermpi.err("Error in command line: {}\n", cmd_args.errorMessage());
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  else if (showhelp) 
    {
      parser.writeToStream(std::cout);
      MPI_Finalize();
    }
  else
    {
      std::stringstream ss;
      ss <<
	"outfile        : " << outfile << std::endl <<
	"latticefile    : " << latticefile << std::endl <<
	"couplingfile   : " << couplingfile << std::endl <<
	"ensemble       : " << ensemble << std::endl <<
	"temperaturefile: " << temperaturefile << std::endl <<
	"np             : " << np << std::endl <<
	"mu             : " << mu << std::endl <<
	"seed           : " << seed << std::endl <<
	"mintemperature : " << mintemperature << std::endl <<
	"precision      : " << precision << std::endl <<
	"iters          : " << iters << std::endl <<
	"verbosity      : " << verbosity << std::endl <<
	"-----" << std::endl;
      lila::loggermpi.out(1, ss.str());      
    }
}

#endif
