#ifndef HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDTHERMO_OPTIONS_H_
#include <string>
#include "clara.hpp"

#include <string>
#include "clara.hpp"

void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, std::string& method, std::string& temperaturefile, int& seed, double& precision, int& neval, int& iters, int& nup, int& ndown, int& np, bool& writeevals, int& argc, char** argv)
{

  bool showhelp = false;
  auto parser =
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(couplingfile, "couplingfile")["-c"]["--couplingfile"]("name of couplingfile (default: latticefile)") |
    clara::Opt(method, "method")["-m"]["--method"]("method to use, either exact or TPQ (default: exact)") |
    clara::Opt(temperaturefile, "temperaturefile")["-t"]["--temperaturefile"]("name of temperaturefile") |
    clara::Opt(seed, "seed")["-s"]["--seed"]("random seed (optional, default 42)") |
    clara::Opt(precision, "precision")["-p"]["--precision"]("precision of Lanczos procedure (optional, default 1e-12)") |
    clara::Opt(neval, "neval")["-n"]["--neval"]("number of eigenvalue to converge (optional, default 1)") |
    clara::Opt(iters, "iters")["-i"]["--iters"]("maximum number of ground state Lanczos iterations performed (optional, default 1000)") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
    clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
    clara::Opt(np, "np")["-p"]["--np"]("number of particles (optional, default half filling)") |
    clara::Opt(writeevals)["-w"]["--writeevals"]("flag whether eigenvalues should not be written to output") |
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
	"method         : " << method << std::endl <<
	"temperaturefile: " << temperaturefile << std::endl <<
	"seed           : " << seed << std::endl <<
	"precision      : " << precision << std::endl <<
	"neval          : " << neval << std::endl <<
	"iters          : " << iters << std::endl <<
	"nup            : " << nup << std::endl <<
	"ndown          : " << ndown << std::endl <<
	"np             : " << np << std::endl <<
	"writeevals     : " << writeevals << std::endl <<
	"-----" << std::endl;
    }
}

#endif
