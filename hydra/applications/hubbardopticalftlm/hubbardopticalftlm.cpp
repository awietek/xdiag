#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "hubbardopticalftlm.options.h"

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
  int iters = 100;
  int verbosity = 1;
  int seed = 1;

  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, iters, verbosity, seed, argc, argv);


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

  /////////////////////////////////////////////////
  // Create Hamiltonian
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

  // Create infrastructure for Hubbard model
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto model = HubbardModelMPI<complex,uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time init: %3.4f\n", t2-t1); 

  // Define multiplication function
  auto multiply_hamiltonian = [&model, &mpi_rank, &verbosity](const VectorMPI<complex>& v, 
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
  //
  /////////////////////////////////////////////////


  /////////////////////////////////////////////////
  // Create current operator
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
  //
  //////////////////////////////////////////////////


  //////////////////////////////////////////////////
  // Create random start state
  int random_seed = seed + 1234567*mpi_rank;
  uint64 local_dim = model.local_dim();
  VectorMPI<complex> startstate(local_dim);
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, random_seed);
  Random(startstate, gen, true);
  Normalize(startstate);
  //
  //////////////////////////////////////////////////


  //////////////////////////////////////////////////
  // Create random start state
  VectorMPI<complex> currstartstate = startstate;
  current.apply_hamiltonian(startstate, currstartstate);
  //
  //////////////////////////////////////////////////


  //////////////////////////////////////////////////
  // Perform Lanczos on start state and current start state
  auto converged = 
    [iters](const lila::Tmatrix<double>& tmat, double beta) {
    (void)beta;
    return LanczosConvergedFixed(tmat, iters);
  };

  std::vector<Vector<complex>> linear_combinations;
  for (int i=0; i< iters; ++i)
    {
      auto lin = Zeros<complex>(iters);
      lin(i) = 1.;
      linear_combinations.push_back(lin);
    }

  auto v0 = startstate;
  auto start_res = Lanczos(multiply_hamiltonian, v0, converged, linear_combinations);
  auto alphas_psi = start_res.tmatrix.diag();
  auto betas_psi = start_res.tmatrix.offdiag();
  auto eigenvalues_psi = start_res.eigenvalues;
  auto& psis = start_res.vectors; 

  v0 = currstartstate;
  Normalize(v0);
  auto current_start_res = Lanczos(multiply_hamiltonian, v0, converged, linear_combinations);
  auto alphas_psi_tilde = current_start_res.tmatrix.diag();
  auto betas_psi_tilde = current_start_res.tmatrix.offdiag();
  auto eigenvalues_psi_tilde = current_start_res.eigenvalues;
  auto& psis_tilde = current_start_res.vectors;
  
  // Compute overlap with start vectors
  Vector<complex> psis_dot_r(iters);
  Vector<complex> psis_tilde_dot_r(iters);
  for (int i=0; i<iters; ++i)
    {
      psis_dot_r(i) = Dot(startstate, psis[i]);
      psis_tilde_dot_r(i) = Dot(psis_tilde[i], currstartstate);
      if (mpi_rank == 0)
	{
	  printf("psis_dot_r(%d)= %.17g %.17g\n", i, std::real(psis_dot_r(i)), std::imag(psis_dot_r(i)));
	  printf("psis_tilde_dot_r(%d)= %.17g %.17g\n\n", i, std::real(psis_tilde_dot_r(i)), std::imag(psis_tilde_dot_r(i)));
	}
    }


  // Compute overlap with start vectors
  Matrix<complex> psis_A_psis_tilde(iters, iters);
  for (int i=0; i<iters; ++i)
    {
      auto tmp = psis[i];
      current.apply_hamiltonian(psis[i], tmp);
      if (mpi_rank == 0)
	printf("computing matrix line %d\n", i); 
	
      for (int j=0; j<iters; ++j)
	  psis_A_psis_tilde(i,j) = Dot(tmp, psis_tilde[j]);
    }

   if (outfile != "")
    {
      if (mpi_rank == 0)
	{
	  std::stringstream line;
	  line << std::setprecision(20);
	  line << "# alphas_psi betas_psi eigenvalues_psi " 
	       << "alphas_psi_tilde betas_psi_tilde eigenvalues_psi_tilde " 
	       << "psis_dot_r psis_tilde_dot_r\n";
	  betas_psi.push_back(0.);
	  betas_psi_tilde.push_back(0.);
	  for (int i=0; i<iters; ++i)
	    line << alphas_psi(i) << " " << betas_psi(i) << " " << eigenvalues_psi(i) << " "
		 << alphas_psi_tilde(i) << " " << betas_psi_tilde(i) << " " << eigenvalues_psi_tilde(i) << " "
		 << std::real(psis_dot_r(i)) << "+" << std::imag(psis_dot_r(i)) << "j " 
		 << std::real(psis_tilde_dot_r(i)) << "+" << std::imag(psis_tilde_dot_r(i))<< "j\n";
	  of << line.str();
	  line.str("");

	  // Write current matrix
	  std::ofstream ofmat;
	  ofmat.open(outfile + std::string(".matrix"), std::ios::out);
	  if(ofmat.fail()) 
	    {
	      if (mpi_rank == 0)
		{
		  std::cerr << "Error in opening matrix outfile: " 
			    << "Could not open file with filename ["
			    << outfile << "] given. Abort." << std::endl;
		  MPI_Abort(MPI_COMM_WORLD, 1);
		}
	    }
	  else
	    {
	      std::stringstream line;
	      line << std::setprecision(20);
	      for (auto i : psis_A_psis_tilde.rows())
		{
		  for (auto j : psis_A_psis_tilde.cols())
		    line << std::real(psis_A_psis_tilde(i,j)) << "+"
			 << std::imag(psis_A_psis_tilde(i,j)) << "j ";
		  line << "\n";
		}
	      ofmat << line.str();
	      line.str("");
	    }

	}

    }




  MPI_Finalize();
  return EXIT_SUCCESS;
}
