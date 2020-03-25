#include "catch.hpp"

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

using namespace hydra::all;
using namespace lila;

void test_tjmodelmpi(hydra::operators::BondList bondlist, 
			  hydra::operators::Couplings couplings, 
			  hydra::hilbertspaces::hubbard_qn qn, 
			  double e0)
{
  auto H = TJModelMPI<double>(bondlist, couplings, qn);
  auto multiply_H = 
    [&H](const VectorMPI<double>& v, VectorMPI<double>& w) 
    { H.apply_hamiltonian(v, w); };
  // Create normal distributed random start state
  VectorMPI<double> startstate(H.local_dim());
  normal_dist_t<double> dist(0., 1.);
  normal_gen_t<double> gen(dist, 42);
  Random(startstate, gen, true);
  Normalize(startstate);  

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply_H, startstate, 1e-12,
				0, "Ritz");
  // printf("e0real: %f, e0comp: %f, diff: %e\n", e0, res.eigenvalues(0), e0 - res.eigenvalues(0));
  REQUIRE(std::abs(e0 - res.eigenvalues(0)) < 1e-6);
}

void test_tjmodelmpi_complex(hydra::operators::BondList bondlist, 
			     hydra::operators::Couplings couplings, 
			     hydra::hilbertspaces::hubbard_qn qn)
{
  auto model = TJModel<complex>(bondlist, couplings, qn);
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  auto full_eigs = lila::EigenvaluesSym(H);


  auto HMPI = TJModelMPI<complex>(bondlist, couplings, qn);
  auto multiply_H = 
    [&HMPI](const VectorMPI<complex>& v, VectorMPI<complex>& w) 
    { HMPI.apply_hamiltonian(v, w); };
  // REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);

  // Create normal distributed random start state
  VectorMPI<complex> startstate(HMPI.local_dim());
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, 42);
  Random(startstate, gen, true);
  Normalize(startstate);  

  // Run Lanczos
  auto res = LanczosEigenvalues(multiply_H, startstate, 1e-12,
				0, "Ritz");
  // std::cout << full_eigs(0) << " " << res.eigenvalues(0) << "\n";
  REQUIRE(std::abs(full_eigs(0) - res.eigenvalues(0)) < 1e-6);
}


TEST_CASE( "TJModelMPI test", "[TJModelMPI]" ) {


  // six site tJ model
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 6;
    for (int s=0; s<n_sites; ++s)
      {
	bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
	bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
      }
   
    // Test t=1.0, J=1.0
    {
      couplings["T"] = 1.0;
      couplings["J"] = 1.0;
      std::vector<double> e0s;
      std::vector<hubbard_qn> qns; 
     
      qns.push_back({0, 0}); e0s.push_back(0.000000000);
      qns.push_back({0, 1}); e0s.push_back(-1.99999999);
      qns.push_back({0, 2}); e0s.push_back(-2.960813110);
      qns.push_back({0, 3}); e0s.push_back(-3.79610527);
      qns.push_back({0, 4}); e0s.push_back(-2.46081311);
      qns.push_back({0, 5}); e0s.push_back(-0.99999999);
      qns.push_back({0, 6}); e0s.push_back(1.500000000);
      qns.push_back({1, 1}); e0s.push_back(-3.61222054);
      qns.push_back({1, 2}); e0s.push_back(-4.04537829);
      qns.push_back({1, 3}); e0s.push_back(-4.10768318);
      qns.push_back({1, 4}); e0s.push_back(-2.42705097);
      qns.push_back({1, 5}); e0s.push_back(-0.49999999);
      qns.push_back({2, 2}); e0s.push_back(-4.16447847);
      qns.push_back({2, 3}); e0s.push_back(-3.52922048);
      qns.push_back({2, 4}); e0s.push_back(-2.11803398);
      qns.push_back({3, 3}); e0s.push_back(-2.80277563);

      for (int i=0; i<(int)qns.size(); ++i)
	test_tjmodelmpi(bondlist, couplings, qns[i], e0s[i]);
    }

    // Test t=1.0, J=0.0
    {
      couplings["T"] = 1.0;
      couplings["J"] = 0.0;
      std::vector<double> e0s;
      std::vector<hubbard_qn> qns; 

      qns.push_back({0, 0}); e0s.push_back(0.000000000);
      qns.push_back({0, 1}); e0s.push_back(-1.99999999);
      qns.push_back({0, 2}); e0s.push_back(-3.00000000);
      qns.push_back({0, 3}); e0s.push_back(-4.00000000);
      qns.push_back({0, 4}); e0s.push_back(-2.99999999);
      qns.push_back({0, 5}); e0s.push_back(-2.00000000);
      qns.push_back({0, 6}); e0s.push_back(0.000000000);
      qns.push_back({1, 1}); e0s.push_back(-3.46410161);
      qns.push_back({1, 2}); e0s.push_back(-3.99999999);
      qns.push_back({1, 3}); e0s.push_back(-3.46410161);
      qns.push_back({1, 4}); e0s.push_back(-1.99999999);
      qns.push_back({1, 5}); e0s.push_back(0.000000000);
      qns.push_back({2, 2}); e0s.push_back(-3.46410161);
      qns.push_back({2, 3}); e0s.push_back(-1.99999999);
      qns.push_back({2, 4}); e0s.push_back(0.000000000);
      qns.push_back({3, 3}); e0s.push_back(0.000000000);

      for (int i=0; i<(int)qns.size(); ++i)
	test_tjmodelmpi(bondlist, couplings, qns[i], e0s[i]);
    }

  }


  // test if complex tJ model MPI gives same e0 as fullED
  {  
    Couplings couplings;
    couplings["T"] = complex(1.0, 1.0);
    couplings["J"] = 1.0;

    // Chains of length 3,4,5,6
    std::vector<int> Ls = {3, 4, 5, 6};
    for (auto L : Ls)
      {
	BondList bondlist;
	for (int s=0; s<L; ++s)
	  {
	    bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % L});
	    bondlist << Bond("HEISENBERG", "J", {s, (s+1) % L});
	  }

	for (int nup=0; nup<=L; ++nup)
	  for (int ndn=0; ndn<L - nup; ++ndn)
	    {
	      hubbard_qn qn = {nup, ndn};
	      test_tjmodelmpi_complex(bondlist, couplings, qn);
	    }
      }

    // Square 3x3
    {
      BondList bondlist;
      bondlist << Bond("HUBBARDHOP", "T", {0, 1});
      bondlist << Bond("HUBBARDHOP", "T", {1, 2});
      bondlist << Bond("HUBBARDHOP", "T", {2, 0});
      bondlist << Bond("HUBBARDHOP", "T", {3, 4});
      bondlist << Bond("HUBBARDHOP", "T", {4, 5});
      bondlist << Bond("HUBBARDHOP", "T", {5, 3});
      bondlist << Bond("HUBBARDHOP", "T", {6, 7});
      bondlist << Bond("HUBBARDHOP", "T", {7, 8});
      bondlist << Bond("HUBBARDHOP", "T", {8, 6});
      bondlist << Bond("HUBBARDHOP", "T", {0, 3});
      bondlist << Bond("HUBBARDHOP", "T", {3, 6});
      bondlist << Bond("HUBBARDHOP", "T", {6, 0});
      bondlist << Bond("HUBBARDHOP", "T", {1, 4});
      bondlist << Bond("HUBBARDHOP", "T", {4, 7});
      bondlist << Bond("HUBBARDHOP", "T", {7, 1});
      bondlist << Bond("HUBBARDHOP", "T", {2, 5});
      bondlist << Bond("HUBBARDHOP", "T", {5, 8});
      bondlist << Bond("HUBBARDHOP", "T", {8, 2});
      bondlist << Bond("HEISENBERG", "J", {0, 1});
      bondlist << Bond("HEISENBERG", "J", {1, 2});
      bondlist << Bond("HEISENBERG", "J", {2, 0});
      bondlist << Bond("HEISENBERG", "J", {3, 4});
      bondlist << Bond("HEISENBERG", "J", {4, 5});
      bondlist << Bond("HEISENBERG", "J", {5, 3});
      bondlist << Bond("HEISENBERG", "J", {6, 7});
      bondlist << Bond("HEISENBERG", "J", {7, 8});
      bondlist << Bond("HEISENBERG", "J", {8, 6});
      bondlist << Bond("HEISENBERG", "J", {0, 3});
      bondlist << Bond("HEISENBERG", "J", {3, 6});
      bondlist << Bond("HEISENBERG", "J", {6, 0});
      bondlist << Bond("HEISENBERG", "J", {1, 4});
      bondlist << Bond("HEISENBERG", "J", {4, 7});
      bondlist << Bond("HEISENBERG", "J", {7, 1});
      bondlist << Bond("HEISENBERG", "J", {2, 5});
      bondlist << Bond("HEISENBERG", "J", {5, 8});
      bondlist << Bond("HEISENBERG", "J", {8, 2});
      
      int n_sites = bondlist.n_sites();
      for (int nup=0; nup<=n_sites; ++nup)
	for (int ndn=0; ndn<n_sites - nup; ++ndn)
	  {
	    hubbard_qn qn = {nup, ndn};
	    test_tjmodelmpi_complex(bondlist, couplings, qn);
	  }
    }

  }

}
