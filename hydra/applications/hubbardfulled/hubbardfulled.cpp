#include <cstdlib>

#include <mpi.h>

#include <lila/all.h>
#include <hydra/all.h>
#include <lime/all.h>

lila::Logger lg;

#include "hubbardfulled.options.h"

template <class coeff_t>
void run_real_complex(std::string real_complex,
		      hydra::all::BondList bondlist,
		      hydra::all::Couplings couplings,
		      hydra::all::hubbard_qn qn,
          bool fulldiag,
		      std::string outfile)
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;
  
  auto file = lime::FileH5(outfile, "w!");
  
  lg.out(1, "Creating {} Hubbard matrixfor n_upspins={}, n_downspins={}...\n",
	 real_complex, qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto model = HubbardModel<coeff_t>(bondlist, couplings, qn);
  auto H = model.matrix();
  double t2 = MPI_Wtime();
  lg.out(1, "done. time: {} secs\n", t2-t1); 
  lg.out(1, "dim: {}\n", FormatWithCommas(model.dim())); 
  
  lg.out(1, "Diagonalizing...\n",
	     qn.n_upspins, qn.n_downspins);
  if (fulldiag == 0) {
    t1 = MPI_Wtime();
    auto eigenvalues = EigenvaluesSym(H);
    t2 = MPI_Wtime();
    lg.out(1, "done. time: {} secs\n", t2-t1); 

    file["Eigenvalues"] = eigenvalues;
    file["Dimension"] = model.dim();
    file.close();

    lg.out(1, "E0: {}\n", eigenvalues(0));
  
  } else {
    t1 = MPI_Wtime();
    auto diagResults = EigenSym(H);
    t2 = MPI_Wtime();
    lg.out(1, "done. time: {} secs\n", t2-t1); 

    file["Eigenvalues"] = diagResults.eigenvalues;
    file["Eigenvectors"] = diagResults.eigenvectors;
    file["Dimension"] = model.dim();
    file.close();

    lg.out(1, "E0: {}\n", diagResults.eigenvalues(0));
  }
}

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  int nup = -1;
  int ndown = -1;
  bool fulldiag = 0;
  int verbosity = 1;

  parse_cmdline(outfile, latticefile, couplingfile, nup, ndown, fulldiag, verbosity, argc, argv);
  lg.set_verbosity(verbosity);  
  
  check_if_files_exists({latticefile, couplingfile});

  // Create Hamiltonian
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);

  // Create infrastructure for Hubbard model
  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};

  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist, couplings, qn, fulldiag, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist, couplings, qn, fulldiag, outfile);
  
  return EXIT_SUCCESS;
}
