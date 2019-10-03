#ifndef HYDRA_APPLICATIONS_HUBBARDOPTICAL_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDOPTICAL_OPTIONS_H_
#include <string>
#include <mpi.h>
#include "clara.hpp"

void parse_cmdline(std::string& outfile, std::string& latticefile, std::string& couplingfile, std::string& corrfile, int& nup, int& ndown, int& iters, int& verbosity, int& seed, int& argc, char** argv)
  {

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  bool showhelp = false;
  auto parser =
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(couplingfile, "couplingfile")["-t"]["--couplingfile"]("name of couplingfile (default: latticefile)") |
    clara::Opt(corrfile, "corrfile")["-c"]["--corrfile"]("name of file describing the current operator") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
    clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
    clara::Opt(iters, "iters")["-i"]["--iters"]("maximum number of ground state Lanczos iterations performed (optional, default 1000)") |
    clara::Opt(verbosity, "verbosity")["-v"]["--verbosity"]("verbosity level, one of 0, 1 ,2 (optional, default 1)") |
    clara::Opt(seed, "seed")["-s"]["--seed"]("random seed") |
    clara::Help(showhelp);


  auto cmd_args = parser.parse(clara::Args(argc,argv));
  if( !cmd_args ) 
    {
      if (mpi_rank==0) std::cerr << "Error in command line: " << cmd_args.errorMessage() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  else if (showhelp) 
    {
      parser.writeToStream(std::cout);
      MPI_Finalize();
    }
  else
    {
      if ((verbosity >= 1) && (mpi_rank == 0))
	{
	  std::cout <<
	    "outfile     : " << outfile << std::endl <<
	    "latticefile : " << latticefile << std::endl <<
	    "couplingfile: " << couplingfile << std::endl <<
	    "corrfile    : " << corrfile << std::endl <<
	    "nup         : " << nup << std::endl <<
	    "ndown       : " << ndown << std::endl <<
	    "iters       : " << iters << std::endl <<
	    "verbosity   : " << verbosity << std::endl <<
	    "seed        : " << seed << std::endl <<
	    "-----" << std::endl;
	}
    }
}

#endif
