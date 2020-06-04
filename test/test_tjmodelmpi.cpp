#include "catch.hpp"

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "testcases_tjmodel.h"

using namespace hydra::all;
using namespace lila;

// Checks whether tjmodelmpi gives same ground state enery as tjmodel
template <class coeff_t>
void test_tjmodelmpi(BondList bondlist,  Couplings couplings)
{
  int n_sites = bondlist.n_sites();
  
  for (int nup=0; nup<=n_sites; ++nup)
    for (int ndn=0; ndn<=n_sites-nup; ++ndn)
      {
	hubbard_qn qn = {nup, ndn};
	auto model = TJModel<complex>(bondlist, couplings, qn);
	auto HM = model.matrix();
	REQUIRE(lila::close(HM, lila::Herm(HM)));
	auto full_eigs = lila::EigenvaluesSym(HM);
	
	auto H = TJModelMPI<coeff_t>(bondlist, couplings, qn);
	auto multiply_H = 
	  [&H](const VectorMPI<coeff_t>& v, VectorMPI<coeff_t>& w) 
	  { H.apply_hamiltonian(v, w); };
	// Create normal distributed random start state
	VectorMPI<coeff_t> startstate(H.local_dim());
	normal_dist_t<coeff_t> dist(0., 1.);
	normal_gen_t<coeff_t> gen(dist, 42);
	Random(startstate, gen, true);
	Normalize(startstate);  

	// Run Lanczos
	auto res = LanczosEigenvalues(multiply_H, startstate, 1e-12,
				      0, "Ritz");
	// LilaPrint(full_eigs(0));
	// LilaPrint(res.eigenvalues(0));
	REQUIRE(std::abs(full_eigs(0) - res.eigenvalues(0)) < 1e-8);
      }
}


TEST_CASE( "TJModelMPI", "[TJModel]" ) {
  using namespace hydra::tjtestcases;
  BondList bondlist;
  Couplings couplings;

  
  ////////////////////////////////////////////////
  // Cross-checks with various heisenberg models
  printf("TJModelMPI: Heisenberg triangle test, N=3\n");
  std::tie(bondlist, couplings) = heisenberg_triangle();
  test_tjmodelmpi<double>(bondlist, couplings);

  
  ///////////////////////////////////////////////////
  // Test all-to-all random coupling Heisenberg
  for (int n_sites=3; n_sites<8; ++n_sites)
    {
      printf("TJModelMPI: Heisenberg random all-to-all test, N=%d\n", n_sites);
      std::tie(bondlist, couplings) = heisenberg_alltoall(n_sites);
      test_tjmodelmpi<double>(bondlist, couplings);
    }

  
  // ///////////////////////////////////////////////////
  // // Kagome 3hexagons
  // printf("TJModelMPI: Heisenberg kagome test, N=15\n");
  // std::tie(bondlist, couplings) = heisenberg_kagome15();
  // test_tjmodelmpi<double>(bondlist, couplings, 4);

  
  // ///////////////////////////////////////////////////
  // // Kagome 6hexagons(outer)
  // printf("TJModelMPI: Heisenberg kagome test, N=39\n");
  // std::tie(bondlist, couplings) = heisenberg_kagome39();
  // test_tjmodelmpi<double>(bondlist, couplings);


  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModelMPI: free fermion random all-to-all test, N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);
      test_tjmodelmpi<double>(bondlist, couplings);
    }
  

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (complex)
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModelMPI: free fermion random all-to-all test (cplx), N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);
      test_tjmodelmpi<complex>(bondlist, couplings);
    }

  ///////////////////////////////////
  // six site tJ model with t=1, J=1
  printf("TJModelMPI: six-site chain test, t=1.0, J=1.0, N=6\n");
  std::tie(bondlist, couplings) = tJchain(6, 1.0, 1.0);
  test_tjmodelmpi<double>(bondlist, couplings);


  ///////////////////////////////////
  // six site tJ model with t=1, J=0
  printf("TJModelMPI: six-site chain test, t=1.0, J=1.0, N=6\n");
  std::tie(bondlist, couplings) = tJchain(6, 1.0, 0.0);
  test_tjmodelmpi<double>(bondlist, couplings);

  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of chains by comparing to ALPS results
  std::vector<int> Ls = {3, 4, 5, 6};
  for (auto L : Ls)
    {
      printf("TJModelMPI: ALPS full spectrum test, chain N=%d\n", L);
      std::tie(bondlist, couplings) = tJchain(L, 1.0, 1.0);
      test_tjmodelmpi<double>(bondlist, couplings);
    }

  // Square 2x2
  printf("TJModelMPI: ALPS full spectrum test, square 2x2\n");
  std::tie(bondlist, couplings) = square2x2(1.0, 1.0);
  test_tjmodelmpi<double>(bondlist, couplings);


  // Square 3x3
  printf("TJModelMPI: ALPS full spectrum test, square 3x3\n");
  std::tie(bondlist, couplings) = square2x2(1.0, 1.0);
  test_tjmodelmpi<double>(bondlist, couplings);

  // test if complex tJ chain gives Hermitian matrix
  Ls = {3, 4, 5, 6};
  for (auto L : Ls)
    {
      printf("TJModelMPI: complex hermitecity test, chain N=%d\n", L);
      std::tie(bondlist, couplings) = tJchain(L, 1.0, 1.0);
      couplings["T"] = complex(1.0, 1.0);
      test_tjmodelmpi<complex>(bondlist, couplings);
    }

  // test if complex tJ all-to-all gives Hermitian matrix
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModelMPI: complex hermitecity test, all-to-all N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);
      test_tjmodelmpi<complex>(bondlist, couplings);
    }
}
