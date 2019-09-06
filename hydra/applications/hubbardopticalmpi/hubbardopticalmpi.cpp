#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "hubbardopticalmpi.options.h"
#include "hubbardopticaliterationsmpi.h"

#include <iomanip>
#include <locale>

template<class T>
std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModelMPI;
  using namespace hydra::operators;
  using hydra::dynamics::continued_fraction;
  using namespace lila;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);


  // Get input parameters
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string corrfile;
  int nup = -1;
  int ndown = -1;
  double precision = 1e-12;
  int iters = 1000;
  double dynprecision = 1e-12;
  int dyniters = 1000;
  int verbosity = 1;
  int lobpcgbands = 1;
  int seed = 0;
  double temperature = 0;

  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, precision, 
		iters, dynprecision, dyniters, verbosity, lobpcgbands, seed, temperature, 
		argc, argv);

  if ((verbosity >= 1) && (mpi_rank == 0))  
    printf("Using %d MPI tasks\n", mpi_size);


  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  if (mpi_rank == 0)
	    {
	      std::cerr << "Error in opening outfile: " 
			<< "Could not open file with filename ["
			<< outfile << "] given. Abort." << std::endl;
	      MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
    }


  // Parse bondlist and couplings from file
  BondList bondlist = read_bondlist(latticefile);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
  if (couplingfile == "") couplingfile = latticefile;
  Couplings couplings = read_couplings(couplingfile);

  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      for (auto bond : hopping_list)
	printf("hopping %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto bond : interaction_list)
	printf("interaction %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto c : couplings)
	printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), std::imag(c.second));
    }

  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  // Parse the corrfile
  std::vector<std::pair<int, int>> correlation_list;
  if (corrfile == "") 
    {
      if (mpi_rank == 0)
	    {
	      std::cerr << "Error in opening corr: " 
			<< "Could not open file with filename ["
			<< corrfile << "] given. Abort." << std::endl;
	      MPI_Abort(MPI_COMM_WORLD, 3);
	    }
    }  
  BondList curr_bondlist = read_bondlist(latticefile);
  

  // Create infrastructure for Hubbard model
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto model = HubbardModelMPI<complex,uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time init: %3.4f\n", t2-t1); 

  // Create infrastructure for current operator
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Creating current operator for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  t1 = MPI_Wtime();
  auto current = HubbardModelMPI<complex,uint32>(curr_bondlist, couplings, qn);
  t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time curr init: %3.4f\n", t2-t1); 


  // Get the TPQ / ground state
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Starting TPQ / ground state Lanczos procedure ...\n");
  int random_seed = 42 + 1234567*mpi_rank;

  auto multiply = [&model, &mpi_rank, &verbosity](const VectorMPI<complex>& v, 
				      VectorMPI<complex>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    model.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    if ((verbosity >= 2) && (mpi_rank == 0))
      printf("iter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };

  MPI_Barrier(MPI_COMM_WORLD);
  uint64 dim = model.dim();
  if ((verbosity >= 1) && (mpi_rank == 0)) printf("dim %s\n", FormatWithCommas(dim).c_str()); 
  MPI_Barrier(MPI_COMM_WORLD);

  uint64 local_dim = model.local_dim();
  VectorMPI<complex> groundstate(local_dim);
  Vector<double> eigs;
  normal_dist_t<complex> dist(1., 0.);
  normal_gen_t<complex> gen(dist, random_seed);

  // Get the groundstate 
  if (temperature == 0)
    {
      // Lanczos algorithm for ground state (unstable)
      if (lobpcgbands == 0)
	{
	  if ((verbosity >= 1) && (mpi_rank == 0))
	    printf("Using Lanczos algorithm for ground state calculation\n");
	  
	  auto res = LanczosEigenvectors(multiply, groundstate, gen, true, precision, {0}, "Ritz");
	  groundstate = res.vectors[0];
	}

      // LOBPCG algorithm for ground state (prefered)
      else
	{
	  if ((verbosity >= 1) && (mpi_rank == 0))
	    {
	      printf("Using LOBPCG algorithm for ground state calculation\n");
	      printf("Using %d bands\n", lobpcgbands);
	    }
	  std::vector<VectorMPI<complex>> vs;
	  for (int i=0; i<lobpcgbands; ++i)
	    {
	      VectorMPI<complex> v(local_dim);
	      Random(v, gen);
	      vs.push_back(v);
	    }
	  auto res = Lobpcg(multiply, vs, precision, iters);
	  eigs = res.eigenvalues;
	  groundstate = res.eigenvectors[0];
	}
    }
  // Get the TPQ state
  else 
    {
      double beta = 1. / temperature;
      Random(groundstate, gen);
      ExpSymVInplace(multiply, groundstate, -beta, precision);  // TODO: get gse
    }

  double norm = lila::real(Dot(groundstate, groundstate));
  groundstate /= (complex)sqrt(norm);
  if ((verbosity >= 1) && (mpi_rank == 0)) printf("Done\n");

  // Debug start
  auto w = groundstate;
  multiply(groundstate, w);
  double e0 = lila::real(Dot(w, groundstate));
  if (mpi_rank == 0)
    {
      printf("eig0: %.17g\n", eigs(0));
      printf("norm: %.17g, e0: %.17g\n", norm, e0);
      // if (!(std::abs(norm - 1) < 1e-12))
      // 	{
      // 	  for (auto a : lzs.alphas())
      // 	    printf("%.17g, ", a);
      // 	  printf("\n");
      // 	  for (auto a : lzs.betas())
      // 	    printf("%.17g, ", a);
      // 	  printf("\n");
      // 	}
    }
  // Debug end

  
  // Write gs energy to outfile
  if (outfile != "")
    {
      if (mpi_rank == 0)
	{
	  std::stringstream line;
	  line << std::setprecision(20);
	  line << "# gs energy: " << eigs(0) << "\n";
	  of << line.str();
	  line.str("");
	}
    }

  
  // for (auto corr : correlation_list)

  //   int s1 = corr.first;
  // int s2 = corr.second;
  // auto res = hydra::hubbard_dynamical_iterations_lanczos_mpi
  //   (model, groundstate, s1, s2, ftype, dyniters, dynprecision, verbosity);

  // // Write to outfile
  // if (mpi_rank == 0)
  //   {
  //     std::stringstream line;

  //     line << "# weight: " << std::setprecision(20) 
  // 	   << res.dyn_weight << "\n";
  //     line << "# alphas betas\n";
  //     of << line.str();
  //     line.str("");

  //     auto alphas = res.alphas;
  //     auto betas = res.betas;
  //     assert(alphas.size() == betas.size());
  //     for (int idx = 0; idx < alphas.size(); ++idx)
  // 	{
  // 	  line << std::setprecision(20)
  // 	       << alphas(idx) << " " << betas(idx) << "\n";
  // 	  of << line.str();
  // 	  line.str("");
  // 	}

  //   }
// }  // algorithm lanczos

MPI_Finalize();
  return EXIT_SUCCESS;
}
