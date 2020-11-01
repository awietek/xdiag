#ifndef HYDRA_APPLICATIONS_HUBBARDFULLONEBODYDENSITY_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDFULLONEBODYDENSITY_OPTIONS_H_
#include <string>
#include <mpi.h>
#include "clara.hpp"

    void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, int& nup, int& ndown, int& verbosity, int& seed, int& argc, char** argv)
  {
  
bool showhelp = false;
 auto parser =
clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
clara::Opt(couplingfile, "couplingfile")["-t"]["--couplingfile"]("name of couplingfile") |
clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
clara::Opt(verbosity, "verbosity")["-v"]["--verbosity"]("verbosity level, one of 0, 1 ,2 (optional, default 1)") |
clara::Opt(seed, "seed")["-s"]["--seed"]("random seed") |
clara::Help(showhelp);


    auto cmd_args = parser.parse(clara::Args(argc,argv));
    if( !cmd_args ) 
      {
         std::cerr << "Error in command line: " << cmd_args.errorMessage() << std::endl;
      }
    else if (showhelp) 
      {
        parser.writeToStream(std::cout);
	MPI_Finalize();
        exit(EXIT_SUCCESS);
      }
    else
      {
	if (verbosity >= 1) 
	  {
	    std::cout <<
	      "outfile     : " << outfile << std::endl <<
	      "latticefile : " << latticefile << std::endl <<
	      "couplingfile: " << couplingfile << std::endl <<
	      "nup         : " << nup << std::endl <<
	      "ndown       : " << ndown << std::endl <<
	      "verbosity   : " << verbosity << std::endl <<
	      "seed        : " << seed << std::endl <<
	      "-----" << std::endl;
	  }
  }
}
#endif
