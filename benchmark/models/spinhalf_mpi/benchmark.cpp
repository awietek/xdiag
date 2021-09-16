#include <hydra/allmpi.h>

#include <unistd.h>

int main(int argc, char *argv[]) {
  using namespace hydra;
  MPI_Init(&argc,&argv);
  LogMPI.set_verbosity(2);

  // int n_sites = 20;
  // std::string lfile = std::string("square.") + std::to_string(n_sites) +
  //                     std::string(".J1J2.fsl.pbc.lat");
  // auto irrep = read_represenation(lfile, "Gamma.C4.A");

  int n_sites = 36;
  std::string lfile = std::string("square.") + std::to_string(n_sites) +
                      std::string(".J1J2.fsl.pbc.lat");

  int n_up = n_sites / 2;
  auto bondlist = read_bondlist(lfile);
  Couplings cpls;
  cpls["J1"] = 1.0;
  tic_mpi();
  auto block = SpinhalfMPI(n_sites, n_up);
  toc_mpi("build block");


  
  tic_mpi();
  double e0 = E0Real(bondlist, cpls, block);
  toc_mpi("Lanczos time");
  LogMPI.out("e0: {}", e0);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
