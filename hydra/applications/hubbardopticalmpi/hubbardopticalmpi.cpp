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

  // Parse input
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
  int seed = 0;
  std::string method = "canonical";
  int lobpcgbands = 1;
  double temperature = 1;

  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown,
		precision, iters, dynprecision, dyniters, verbosity, seed,
		method, lobpcgbands, temperature, argc, argv);

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
  BondList curr_bondlist = read_bondlist(corrfile);
  

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
  for (auto bond : curr_bondlist)
    if ((verbosity >= 1) && (mpi_rank == 0))
      printf("curr %s %d %d\n", bond.coupling().c_str(), 
	     bond.sites()[0], bond.sites()[1]);

  Couplings curr_couplings;
  curr_couplings["C"] = 1;
  auto current = HubbardModelMPI<complex,uint32>(curr_bondlist, curr_couplings, qn);
  t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time curr init: %3.4f\n", t2-t1); 


  // Get the TPQ / ground state
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Starting TPQ / ground state Lanczos procedure ...\n");
  int random_seed = seed + 1234567*mpi_rank;

  // Define multiplication function
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
  VectorMPI<complex> tpqstate(local_dim);
  Vector<double> eigs;
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, random_seed);
  
  // Get the tpqstate
  double e0 = 0;
  double measured_temperature = 0;
  if (method == "groundstate")
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
      tpqstate = res.eigenvectors[0];
      e0 = eigs(0);
      measured_temperature = 0;
    }
  // Get the canonical TPQ state
  else if (method == "canonical") 
    {
      if ((verbosity >= 1) && (mpi_rank == 0))
	printf("Using Lanczos algorithm for canonical TPQ state\n");

      double beta = 1. / temperature;
      Random(tpqstate, gen);
      e0 = ExpSymVInplace(multiply, tpqstate, -beta / 2., precision, true);
      measured_temperature = temperature;
    }
  // Get the microcanonical TPQ state
  else if (method == "microcanonical") 
    {
      if ((verbosity >= 1) && (mpi_rank == 0))
	printf("Using Lanczos algorithm for microcanonical TPQ state\n");

      // First run for Tmatrix
      Random(tpqstate, gen, false);  // do not alter the generator
      int n_eigenvalue = 0;
      auto res_first = 
	LanczosEigenvalues(multiply, tpqstate, precision, n_eigenvalue, "Ritz");
      // Define fixed iterations convergence criterion
      int n_iterations = res_first.tmatrix.size();
      auto converged = 
	[n_iterations](const lila::Tmatrix<double>& tmat, double beta) {
	  (void)beta;
	  return LanczosConvergedFixed(tmat, n_iterations);
	};
    
      // Prepare linear combinations
      int n_eigenvector = (int)temperature;
      auto evecs = Eigen(res_first.tmatrix).eigenvectors;
      std::vector<Vector<complex>> linear_combinations;
      Vector<complex> evc(n_iterations);
      for (int i=0; i<n_iterations; ++i)
	evc(i) = (complex)evecs(i, n_eigenvector);
      linear_combinations.push_back(evc);
      // reset initial vector and rerun with linear combination
      Random(tpqstate, gen, true);
      auto res = Lanczos(multiply, tpqstate, converged, linear_combinations);
      eigs = res.eigenvalues;
      tpqstate = res.vectors[0];
      e0 = eigs(0);
      measured_temperature = 1.23;
    }

  double norm = lila::real(Dot(tpqstate, tpqstate));
  tpqstate /= (complex)sqrt(norm);
  if ((verbosity >= 1) && (mpi_rank == 0)) printf("Done\n");

  auto w = tpqstate;
  multiply(tpqstate, w);
  double energy = lila::real(Dot(w, tpqstate));
  double energy2 = lila::real(Dot(w, w));
  double delta_e = energy2 - energy*energy;
  auto startstate = tpqstate;
  current.apply_hamiltonian(tpqstate, startstate);
  double weight = lila::Norm(startstate);

  if (mpi_rank == 0)
    {
      printf("norm: %.17g, e0: %.17g, e: %.17g, delta_e: %.17g, weight: %.17g\n",
	     norm, e0, energy, delta_e, weight);
    }
  
  // Write energies and weight to outfile
  if (outfile != "")
    {
      if (mpi_rank == 0)
	{
	  std::stringstream line;
	  line << std::setprecision(20);
	  line << "# weight: " << weight << "\n";
	  line << "# energy: " << energy << "\n";
	  line << "# temperature: " << measured_temperature << "\n";
	  line << "# gs energy: " << e0 << "\n";
	  of << line.str();
	  line.str("");
	}
    }
  
  int num_eigenvalue = 1;
  auto lczs_res = LanczosEigenvalues(multiply, startstate, dynprecision, 
				     num_eigenvalue, "Eigenvalues");

  auto eigenvalues = lczs_res.eigenvalues;
  auto alphas = lczs_res.tmatrix.diag();
  auto betas = lczs_res.tmatrix.offdiag();
  if (mpi_rank==0)
    {
      LilaPrint(eigenvalues);
      LilaPrint(alphas);
      LilaPrint(betas);
    }

   if (outfile != "")
    {
      if (mpi_rank == 0)
	{
	  std::stringstream line;
	  line << std::setprecision(20);
	  line << "# alphas betas eigenvalues\n";
	  betas.push_back(0.);
	  for (int i=0; i<alphas.size(); ++i)
	    line << alphas(i) << " " << betas(i) << " " << eigenvalues(i) << "\n";
	  of << line.str();
	  line.str("");
	}
    }
   
  MPI_Finalize();
  return EXIT_SUCCESS;
}
