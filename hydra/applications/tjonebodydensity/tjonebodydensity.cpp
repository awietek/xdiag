#include <cstdlib>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>
#include <lime/all.h>

#include <stdio.h>

lila::LoggerMPI lg;

#include "tjonebodydensity.options.h"

template <class coeff_t>
void run_real_complex(std::string real_complex,
    hydra::all::BondList bondlist,
    hydra::all::Couplings couplings,
    hydra::all::hubbard_qn qn,
    double precision, 
    int iters, std::string algorithm,
     int lobpcgbands,
    int verbosity, int seed,
    std::string outfile)
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  FileH5 file;
  if (mpi_rank == 0) file = lime::FileH5(outfile, "w!");

  lg.out(1, "Creating {} t-J model for n_upspins={}, n_downspins={}...\n",
      real_complex, qn.n_upspins, qn.n_downspins);
  lg.out(1, "Using {} MPI tasks\n", mpi_size);
  
  // Create tJ Hamiltonian 

  int n_sites = bondlist.n_sites();
  double t1 = MPI_Wtime();
  auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  lg.out(1, "time init: {} secs\n", t2-t1); 


  // Define Hamilton multiplication function
  int iter = 0;
  auto multiply_H = 
    [&H, &iter, verbosity](const VectorMPI<coeff_t>& v, VectorMPI<coeff_t>& w) 
    {
      bool verbose = (iter==0) && (verbosity > 0);
      double t1 = MPI_Wtime();
      H.apply_hamiltonian(v, w, verbose);
      double t2 = MPI_Wtime();
      lg.out(2, "iter: {}, time MVM: {}\n", iter, t2-t1); 
      ++iter;
    };


  ///////////////////////////////
  // Get the ground state
  auto local_dim = H.local_dim();
  VectorMPI<coeff_t> groundstate(local_dim);
  Vector<double> eigs;

  int random_seed = seed + 1234567*mpi_rank;
  normal_dist_t<coeff_t> dist(0., 1.);
  normal_gen_t<coeff_t> gen(dist, random_seed);

  // Lanczos algorithm for ground state (unstable)
  if (lobpcgbands == 0)
    {
      lg.out(1, "Starting ground state Lanczos...\n");
      auto res = LanczosEigenvectors(multiply_H, groundstate, gen, true, precision, {0}, "Ritz");
      eigs = res.eigenvalues;
      groundstate = res.vectors[0];
    }

  // LOBPCG algorithm for ground state (prefered)
  else
    {
      lg.out(1, "Starting ground state LOBPCG algorithm...\n");
      lg.out(1, "Using {} bands\n", lobpcgbands);
      
      // Create random start vectors
      std::vector<VectorMPI<coeff_t>> vs;
      for (int i=0; i<lobpcgbands; ++i)
	{
	  VectorMPI<coeff_t> v(local_dim);
	  Random(v, gen);
	  vs.push_back(v);
	}

      // Run LOBPCG
      auto res = Lobpcg(multiply_H, vs, precision, iters);
      eigs = res.eigenvalues;
      groundstate = res.eigenvectors[0];
    }

  coeff_t norm = Dot(groundstate, groundstate);
  groundstate /= sqrt(norm);
  lg.out(1, "Done\n");
  lg.out(1, "Ground state energy: {}\n", eigs(0));
  if (mpi_rank == 0)
  {
    file["GroundStateEnergy"] = eigs(0);
  }

  // Generate one-body density matrix

  VectorMPI<coeff_t> perturbedGroundstate(local_dim);
  lila::Matrix<complex> singleParticleCorrelations;
  singleParticleCorrelations.resize(n_sites, n_sites);

  // First, generate a bondlist with all-to-all coupling
  
  std::vector<Bond> bonds;
  for (int i=0;i<n_sites;i++) {
    for (int j=0;j<n_sites;j++) {
      if (i != j)
      {
        std::vector<int> sites{i, j};
        std::string hopping_label = std::to_string(i) + std::to_string(j);
        bonds.push_back(Bond("HUBBARDHOP", hopping_label, sites));
      } else {
        std::vector<int> sites{i};
        std::string hopping_label = std::to_string(i);
        bonds.push_back(Bond("HUBBARDMU", hopping_label, sites));
      }
    }
  }
  BondList densityBondlist(bonds);
  
  // Measure diagonal elements first

  for (int i=0;i<n_sites;i++) {
    std::map<std::string, complex> coupling_map;
    std::string hopping_label = std::to_string(i);
    coupling_map[hopping_label] = -1;
    Couplings densityCouplings(coupling_map);
    auto hoppingH = TJModelMPI<coeff_t>(densityBondlist, densityCouplings, qn);
    hoppingH.apply_hamiltonian(groundstate, perturbedGroundstate, false);
    singleParticleCorrelations(i, i) = Dot(perturbedGroundstate, groundstate);
  }

  // Then measure off-diagonal elements.
  // If hoppings are complex, it is necessary for the one-body density matrix
  // to measure <c_i^\dagger c_j> and <c_j^\dagger c_i> separately.
  // Note: the technical reason this exception is necessary is because complex
  // hoppings won't work if coeff_t is real, since the ground state will be real.
  // This would not be necessary if, in the future,
  // complex Hamiltonians can act on real ground states.

  if (real_complex == std::string("COMPLEX")) {
    for (int i=0;i<n_sites;i++) {
      for (int j=(i+1);j<n_sites;j++) {

        // First measure <c_i^\dagger c_j + c_j^\dagger c_j>
        
        std::map<std::string, complex> coupling_map;
        std::string hopping_label = std::to_string(i) + std::to_string(j);
        coupling_map[hopping_label] = 1;
        Couplings densityCouplings(coupling_map);
        auto hopping_symmetric = TJModelMPI<coeff_t>(densityBondlist, densityCouplings, qn);
        hopping_symmetric.apply_hamiltonian(groundstate, perturbedGroundstate, false);
        complex symmetric_expectation = Dot(perturbedGroundstate, groundstate);

        // Then measure <c_i^\dagger c_j - c_j^\dagger c_i>
        coupling_map[hopping_label] = std::complex<double>(0, 1);
        densityCouplings = Couplings(coupling_map);
        auto hopping_antisymmetric = TJModelMPI<coeff_t>(densityBondlist, densityCouplings, qn);
        hopping_antisymmetric.apply_hamiltonian(groundstate, perturbedGroundstate, false);
        complex antisymmetric_expectation = complex(0, -1)*Dot(perturbedGroundstate, groundstate);
        singleParticleCorrelations(i, j) = 0.5*(symmetric_expectation + antisymmetric_expectation);
        singleParticleCorrelations(j, i) = 0.5*(symmetric_expectation - antisymmetric_expectation);
      }
    }
  } else {
    for (int i=0;i<n_sites;i++) {
      for (int j=(i+1);j<n_sites;j++) {
        std::map<std::string, complex> coupling_map;
        std::string hopping_label = std::to_string(i) + std::to_string(j);
        coupling_map[hopping_label] = 1;
        Couplings densityCouplings(coupling_map);
        auto hoppingH = TJModelMPI<coeff_t>(densityBondlist, densityCouplings, qn);
        hoppingH.apply_hamiltonian(groundstate, perturbedGroundstate, false);
        singleParticleCorrelations(i,j) = 0.5*Dot(perturbedGroundstate, groundstate);
        singleParticleCorrelations(j,i) = 0.5*Dot(perturbedGroundstate, groundstate);
      }
    }
  }

  if (mpi_rank == 0) {
  file["Correlations"] = singleParticleCorrelations;
  file["Hoppings"] = H.single_particle_hopping();
	file.close();
  }
}

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  MPI_Init(&argc, &argv); 

  // Get input parameters
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  int nup = -1;
  int ndown = -1;
  std::string algorithm;
  int iters = 1000;
  double precision = 1e-10;
  int seed = 1;
  int verbosity = 1;
  int lobpcgbands = 1;
  parse_cmdline(outfile, latticefile, couplingfile, nup, ndown, 
		 algorithm, precision, iters, verbosity, lobpcgbands, seed, argc, argv);
  lg.set_verbosity(verbosity);  
  check_if_files_exists({latticefile, couplingfile});
  check_if_contained_in(algorithm, {"lanczos", "bandlanczos"});


  // Parse bondlist, couplings, and corr_bondlist from file
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);


  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist,
        couplings, qn, precision, iters, algorithm,
        lobpcgbands, verbosity, seed, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist,
        couplings, qn, precision, iters, algorithm,
        lobpcgbands, verbosity, seed, outfile);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
  
