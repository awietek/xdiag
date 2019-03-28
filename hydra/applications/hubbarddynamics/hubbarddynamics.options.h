#ifndef HYDRA_APPLICATIONS_HUBBARDDYNAMICS_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDDYNAMICS_OPTIONS_H_
#include <string>
#include "clara.hpp"


void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, std::string& corrfile, int& nup, int& ndown, std::string& fermiontype, std::string& algorithm, double& precision, int& iters, double& dynprecision, int& dyniters, int& verbosity, double& deflationtol, int& argc, char** argv)
  {

bool showhelp = false;
 auto parser =
clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
clara::Opt(couplingfile, "couplingfile")["-t"]["--couplingfile"]("name of couplingfile (default: latticefile)") |
clara::Opt(corrfile, "corrfile")["-c"]["--corrfile"]("name of file containing sites of correlator (optional: default empty, only 0-0 correlator is computed)") |
clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
clara::Opt(fermiontype, "fermiontype")["-f"]["--fermiontype"]("type of fermion (one of cdagdn, cdn) (optinal: default all are computed)") |
clara::Opt(algorithm, "algorithm")["-a"]["--algorithm"]("algorithm used (one of lanczos, bandlanczos) (default: lanczos)") |
clara::Opt(precision, "precision")["-p"]["--precision"]("precision of ground state Lanczos procedure (optional, default 1e-12)") |
clara::Opt(iters, "iters")["-i"]["--iters"]("maximum number of ground state Lanczos iterations performed (optional, default 1000)") |
clara::Opt(dynprecision, "dynprecision")["-q"]["--dynprecision"]("precision of dynamical Lanczos procedure (optional, default 1e-12)") |
clara::Opt(dyniters, "dyniters")["-j"]["--dyniters"]("maximum number of dynamical Lanczos iterations performed (optional, default 1000)") |
clara::Opt(verbosity, "verbosity")["-v"]["--verbosity"]("verbosity level, one of 0, 1 ,2 (optional, default 1)") |
clara::Opt(deflationtol, "deflationtol")["-x"]["--deflationtol"]("tolerance for detecting deflation in Band Lanczos algortihm (optional, default 1e-8)") |
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
      "corrfile    : " << corrfile << std::endl <<
      "nup         : " << nup << std::endl <<
      "ndown       : " << ndown << std::endl <<
      "fermiontype : " << fermiontype << std::endl <<
      "algorithm   : " << algorithm << std::endl <<
      "precision   : " << precision << std::endl <<
      "iters       : " << iters << std::endl <<
      "dynprecision: " << dynprecision << std::endl <<
      "dyniters    : " << dyniters << std::endl <<
      "verbosity   : " << verbosity << std::endl <<
      "deflationtol: " << deflationtol << std::endl <<
      "-----" << std::endl;
  }
}

#endif
