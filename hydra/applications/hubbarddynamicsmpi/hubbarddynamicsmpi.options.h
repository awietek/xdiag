#ifndef HYDRA_APPLICATIONS_HUBBARDDYNAMICS_OPTIONS_H_
#define HYDRA_APPLICATIONS_HUBBARDDYNAMICS_OPTIONS_H_
#include <string>
#include "clara.hpp"

void parse_cmdline(double& t, double& U, std::string& outfile, 
		   std::string& latticefile, int& nup, int& ndown, 
		   std::string& fermiontype, int& dyniters, int& argc, char** argv)
{

  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  bool showhelp = false;
  auto parser =
    clara::Opt(t, "t")["-t"]["--t"]("Hubbard t") |
    clara::Opt(U, "U")["-U"]["--U"]("Hubbard U") |
    clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
    clara::Opt(latticefile, "latticefile")["-l"]["--latticefile"]("name of latticefile") |
    clara::Opt(nup, "nup")["-u"]["--nup"]("number of up electrons (optional, default half filling)") |
    clara::Opt(ndown, "ndown")["-d"]["--ndown"]("number of down electrons (optional, default half filling)") |
    clara::Opt(fermiontype, "fermiontype")["-f"]["--fermiontype"]("type of fermion (one of cdagup, cdag, cdagdn, cdn) (optinal: default all are computed)") |
    clara::Opt(dyniters, "dyniters")["-d"]["--dyniters"]("number of dynamical Lanczos iterations performed (optional, default 200)") |
    clara::Help(showhelp);


  auto cmd_args = parser.parse(clara::Args(argc,argv));
  if( !cmd_args ) 
    {
      if (mpi_rank == 0) std::cerr << "Error in command line: " << cmd_args.errorMessage() << std::endl;
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  else if (showhelp) 
    {
      if (mpi_rank == 0) parser.writeToStream(std::cout);
      MPI_Finalize();
      exit(EXIT_SUCCESS);
    }
  else
    {
      if (mpi_rank == 0)
	{
	  std::cout <<
	    "t          : " << t << std::endl <<
	    "U          : " << U << std::endl <<
	    "outfile    : " << outfile << std::endl <<
	    "latticefile: " << latticefile << std::endl <<
	    "nup        : " << nup << std::endl <<
	    "ndown      : " << ndown << std::endl <<
	    "fermiontype: " << fermiontype << std::endl <<
	    "dyniters   : " << dyniters << std::endl <<
	    "-----" << std::endl;
	}
    }
}
#endif
