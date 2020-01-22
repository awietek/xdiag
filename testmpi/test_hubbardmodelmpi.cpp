#include "catch.hpp"

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

using namespace hydra::all;
using namespace lila;

void test_hubbardmodelmpi(hydra::operators::BondList bondlist, 
			  hydra::operators::Couplings couplings, 
			  hydra::hilbertspaces::hubbard_qn qn, 
			  double e0)
{


  auto H = HubbardModelMPI<double>(bondlist, couplings, qn);
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
  printf("e0real: %f, e0comp: %f, diff: %e\n", e0, res.eigenvalues(0), e0 - res.eigenvalues(0));
  REQUIRE(std::abs(e0 - res.eigenvalues(0)) < 1e-10);
}

TEST_CASE( "Add MPI test", "[AddMPI]" ) {

  //////////////////////////////////////
  // Test free fermion chains


  // Two sites free-fermions
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 2;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
    couplings["T"] = 1.0;
    double e0 = -4.0;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // four sites free-fermions
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 4;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
    couplings["T"] = 1.0;
    double e0 = -4.0;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // six sites free-fermions
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 6;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
    couplings["T"] = 1.0;
    double e0 = -8.0;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // eight sites free-fermions
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 8;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % n_sites});
    couplings["T"] = 1.0;
    double e0 = -9.65685424949;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }




  //////////////////////////////////////////////
  // Test antiferromagnetic heisenberg models
  // Two sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 2;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -1.5;  // two bonds on same sites -> 2*3/4=1.5
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // Three sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 3;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -0.75;  // two bonds on same sites -> 2*3/4=1.5
    hubbard_qn qn = {1, 2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }


  // four sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 4;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -2.0;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // five-sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 5;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -1.8680339887498953;
    hubbard_qn qn = {2, 3};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // six sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 6;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -2.8027756377319957;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }

  // eight sites Heisenberg
  {  
    BondList bondlist;
    Couplings couplings;

    int n_sites = 8;
    for (int s=0; s<n_sites; ++s)
      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % n_sites});
    couplings["J"] = 1.0;
    double e0 = -3.651093408937172;
    hubbard_qn qn = {n_sites/2, n_sites/2};
    test_hubbardmodelmpi(bondlist, couplings, qn, e0);
  }



}
