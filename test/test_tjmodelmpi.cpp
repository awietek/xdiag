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


}
