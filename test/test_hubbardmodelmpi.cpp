#include "catch.hpp"

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "testcases_hubbardmodel.h"

using namespace hydra::all;
using namespace lila;

template <class coeff_t>
void test_hubbardmodelmpi(hydra::operators::BondList bondlist, 
			  hydra::operators::Couplings couplings)
{
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  int n_sites = bondlist.n_sites();
  for (int nup=0; nup <= n_sites; ++nup)
    for (int ndn=0; ndn <= n_sites; ++ndn)
      {
	hubbard_qn qn = {nup, ndn};
	
	auto H = HubbardModelMPI<coeff_t>(bondlist, couplings, qn);
	auto multiply_H = 
	  [&H](const VectorMPI<coeff_t>& v, VectorMPI<coeff_t>& w) 
	  { H.apply_hamiltonian(v, w); };
	// Create normal distributed random start state
	VectorMPI<coeff_t> startstate(H.local_dim());
	normal_dist_t<coeff_t> dist(0., 1.);
	normal_gen_t<coeff_t> gen(dist, 42 + mpi_rank);
	Random(startstate, gen, true);
	Normalize(startstate);  
	
	// Run Lanczos
	auto res = LanczosEigenvalues(multiply_H, startstate, 1e-12,
				      0, "Ritz");
	
	auto model = HubbardModel<coeff_t>(bondlist, couplings, qn);
	auto HM = model.matrix();
	auto eigs = lila::EigenvaluesSym(HM);
	auto e0 = eigs(0);
	
	REQUIRE(std::abs(e0 - res.eigenvalues(0)) < 1e-10);
      }
}

TEST_CASE( "HubbardModelMPI", "[HubbardModelMPI]" )
{
  using namespace hydra::hubbardtestcases;

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  BondList bondlist;
  Couplings couplings;
  
  //////////////////////////////////
  // Test two site exact solution
  bondlist << Bond("HUBBARDHOP", "T", {0, 1});
  couplings["T"] = 1.0;
  for (int i=0; i<20; ++i)
    {
      double U = 1.234 * i;
      couplings["U"] = U;
      if (mpi_rank==0) printf("HubbardModelMPI: two-site exact solution test, U=%f\n", U);
      test_hubbardmodelmpi<double>(bondlist, couplings);
    }

  ////////////////////////////////////////////////
  // Cross-checks with various heisenberg models
  if (mpi_rank==0) printf("HubbardModelMPI: Heisenberg triangle test, N=3\n");
  std::tie(bondlist, couplings) = heisenberg_triangle();
  test_hubbardmodelmpi<double>(bondlist, couplings);
  

  ///////////////////////////////////////////////////
  // Test all-to-all random coupling Heisenberg
  for (int n_sites=3; n_sites<8; ++n_sites)
    {
      if (mpi_rank==0)
	{
	  printf("HubbardModelMPI: Heisenberg random all-to-all test, ");
	  printf("N=%d\n", n_sites);
	}
      std::tie(bondlist, couplings) = heisenberg_alltoall(n_sites);
      test_hubbardmodelmpi<double>(bondlist, couplings);
    }
  
  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (real)
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      if (mpi_rank==0)
	{
	  printf("HubbardModelMPI: freefermion random all-to-all test, ");
	  printf("N=%d\n", n_sites);
	}
      std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);
      test_hubbardmodelmpi<double>(bondlist, couplings);
    }  

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (complex)
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      if (mpi_rank==0)
	{
	  printf("HubbardModelMPI: freefermion random all-to-all test (cplx), ");
	  printf("N=%d\n", n_sites);
	}
      std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);
      test_hubbardmodelmpi<complex>(bondlist, couplings);
    }
 
  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of random all-to-all interactions 
  // by comparing to MATLAB results
  if (mpi_rank==0) printf("HubbardModelMPI: MATLAB full spectrum, chain N=4, no Hubbard U\n");
  std::tie(bondlist, couplings) = randomAlltoAll4NoU();
  test_hubbardmodelmpi<double>(bondlist, couplings);

  if (mpi_rank==0) printf("HubbardModelMPI: MATLAB full spectrum, chain N=4, Hubbard U\n");
  std::tie(bondlist, couplings) = randomAlltoAll4();
  test_hubbardmodelmpi<double>(bondlist, couplings);

}
