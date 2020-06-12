#include "catch.hpp"

#include <fstream>
#include <lila/all.h>
#include <hydra/all.h>
#include <random>

#include "testcases_tjmodel.h"

using namespace hydra::all;
using namespace lila;

// Test by comparing to Heisenberg model
void test_tjmodel_heisenberg(hydra::operators::BondList bondlist, 
			     hydra::operators::Couplings couplings,
			     int max_nup=-1)
{
  int n_sites = bondlist.n_sites();
  for (int nup=0; nup<= (max_nup < 0 ? n_sites : std::min(max_nup, n_sites));
       ++nup)
    {
      int ndn = n_sites - nup;
      hydra::hilbertspaces::hubbard_qn qn = {nup, ndn};

      auto HB = HeisenbergModel();
      auto HBM = HB.matrix(bondlist, couplings, nup);
      
      // Create tJ matrix
      auto TJ = TJModel<double>(bondlist, couplings, qn);
      auto TJM = TJ.matrix();
      
      auto etj = lila::EigenvaluesSym(TJM);
      auto ehb = lila::EigenvaluesSym(HBM);

      REQUIRE(lila::close(HBM, lila::Herm(HBM)));
      REQUIRE(lila::close(TJM, lila::Herm(TJM)));
      
      REQUIRE(etj.size() == ehb.size());
      REQUIRE(lila::close(etj, ehb));

      // printf("tjhb nup: %d, ndn: %d, etj(0): %f, ehb(0) %f\n",
      // 	     nup, ndn, etj(0), ehb(0));
    }
}

// Test by comparing to exactly known solution for given quantumnumber
template <class coeff_t=double>
void test_tjmodel_e0(hydra::operators::BondList bondlist, 
		     hydra::operators::Couplings couplings, 
		     hydra::hilbertspaces::hubbard_qn qn, 
		     double e0)
{
  auto model = TJModel<coeff_t>(bondlist, couplings, qn);
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  auto eigs = lila::EigenvaluesSym(H);
  // printf("tje0 t: %f, J: %f, nup: %d, ndn: %d, eigs(0): %f, e0 %f\n",
  // 	 lila::real(couplings["T"]), lila::real(couplings["J"]),
  // 	 qn.n_upspins, qn.n_downspins, eigs(0), e0);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}


// Test by comparing to full spectrum of alps
void test_tjmodel_fullspectrum(hydra::operators::BondList bondlist, 
		       hydra::operators::Couplings couplings,
           bool ninj_term,
		       std::string filename)
{
  int n_sites = bondlist.n_sites();

  // Compute full spectrum in hydra
  lila::Vector<double> all_eigs;
  for (int nup=0; nup<=n_sites; ++nup)
    for (int ndn=0; ndn<=n_sites - nup; ++ndn)
      {
	hubbard_qn qn = {nup, ndn};
	auto model = TJModel<double>(bondlist, couplings, qn);

	// Run Full ED
	auto H = model.matrix(ninj_term);

	if (nup + ndn == n_sites)
	  REQUIRE(lila::close(H, lila::Herm(H)));
	auto eigs = lila::EigenvaluesSym(H);
	for (auto eig : eigs)
	  all_eigs.push_back(eig);
  
      }
  std::sort(all_eigs.begin(), all_eigs.end());

  // Read alps_eigs
  lila::Vector<double> alps_eigs;
  std::ifstream in(filename.c_str());
  if(!in)
    {
      std::cerr << "test_tjmodel.cpp: Cannot open the File : "
		<< filename <<std::endl;
      exit(EXIT_FAILURE);
    }
  std::string str;
  while (std::getline(in, str))
    {
      if(str.size() > 0)
        alps_eigs.push_back(std::stod(str));
    }
  in.close();
  std::sort(alps_eigs.begin(), alps_eigs.end());
  REQUIRE(all_eigs.size() == alps_eigs.size());  
  for (int i=0; i<all_eigs.size(); ++i)
    REQUIRE(close(all_eigs(i), alps_eigs(i)));
}


// Test if complex matrix is hermitian
void test_tjmodel_hermitian(hydra::operators::BondList bondlist, 
			    hydra::operators::Couplings couplings)
{
  int n_sites = bondlist.n_sites();
  for (int nup=0; nup<=n_sites; ++nup)
    for (int ndn=0; ndn<=n_sites-nup; ++ndn)
      {
 	auto model = TJModel<complex>(bondlist, couplings, {nup, ndn});
	auto H = model.matrix();
	REQUIRE(lila::close(H, lila::Herm(H)));
      }
}

TEST_CASE( "TJModel", "[TJModel]" )
{
  using namespace hydra::tjtestcases;
  BondList bondlist;
  Couplings couplings;

  
  ////////////////////////////////////////////////
  // Cross-checks with various heisenberg models
  printf("TJModel: Heisenberg triangle test, N=3\n");
  std::tie(bondlist, couplings) = heisenberg_triangle();
  test_tjmodel_heisenberg(bondlist, couplings);

  
  ///////////////////////////////////////////////////
  // Test all-to-all random coupling Heisenberg
  for (int n_sites=3; n_sites<8; ++n_sites)
    {
      printf("TJModel: Heisenberg random all-to-all test, N=%d\n", n_sites);
      std::tie(bondlist, couplings) = heisenberg_alltoall(n_sites);
      test_tjmodel_heisenberg(bondlist, couplings);
    }

  
  ///////////////////////////////////////////////////
  // Kagome 3hexagons
  printf("TJModel: Heisenberg kagome test, N=15\n");
  std::tie(bondlist, couplings) = heisenberg_kagome15();
  test_tjmodel_heisenberg(bondlist, couplings, 4);

  
  ///////////////////////////////////////////////////
  // Kagome 6hexagons(outer)
  printf("TJModel: Heisenberg kagome test, N=39\n");
  std::tie(bondlist, couplings) = heisenberg_kagome39();
  test_tjmodel_heisenberg(bondlist, couplings, 2);


  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModel: free fermion random all-to-all test, N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall(n_sites);

      // Create single particle matrix
      auto Hs = lila::Zeros<double>(n_sites, n_sites);
      for (auto bond : bondlist)
	{
	  assert(bond.size() == 2);
	  int s1 = bond.sites(0);
	  int s2 = bond.sites(1);
	  auto name = bond.coupling();
	  Hs(s1, s2) = -lila::real(couplings[name]);
	  Hs(s2, s1) = -lila::real(couplings[name]);
	}
      auto seigs = lila::EigenvaluesSym(Hs);
      for (int nup=0; nup <= n_sites; ++nup)
	  {
	    double e0 = 0;
	    for (int i=0; i<nup; ++i)
	      e0 += seigs(i);
	    test_tjmodel_e0(bondlist, couplings, {nup,0}, e0);
	    test_tjmodel_e0(bondlist, couplings, {0,nup}, e0);
	  }
    }
  

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (complex)
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModel: free fermion random all-to-all test (cplx), N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);

      // Create single particle matrix
      auto Hs = lila::Zeros<complex>(n_sites, n_sites);
      for (auto bond : bondlist)
	{
	  assert(bond.size() == 2);
	  int s1 = bond.sites(0);
	  int s2 = bond.sites(1);
	  auto name = bond.coupling();
	  Hs(s1, s2) = -couplings[name];
	  Hs(s2, s1) = -lila::conj(couplings[name]);
	}
      auto seigs = lila::EigenvaluesSym(Hs);
      
      for (int nup=0; nup <= n_sites; ++nup)
	  {
	    double e0 = 0;
	    for (int i=0; i<nup; ++i)
	      e0 += seigs(i);
	    test_tjmodel_e0<complex>(bondlist, couplings, {nup,0}, e0);
	    test_tjmodel_e0<complex>(bondlist, couplings, {0,nup}, e0);
	  }
    }

  ///////////////////////////////////
  // six site tJ model with t=1, J=1
  printf("TJModel: six-site chain test, t=1.0, J=1.0, N=6\n");
  std::tie(bondlist, couplings) = tJchain(6, 1.0, 1.0);

  // Known solutions
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
    test_tjmodel_e0(bondlist, couplings, qns[i], e0s[i]);


  ///////////////////////////////////
  // six site tJ model with t=1, J=0
  printf("TJModel: six-site chain test, t=1.0, J=1.0, N=6\n");
  std::tie(bondlist, couplings) = tJchain(6, 1.0, 0.0);
  
  // Known solutions
  e0s.clear();
  qns.clear();
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
    test_tjmodel_e0(bondlist, couplings, qns[i], e0s[i]);
   
  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of chains by comparing to ALPS results
  std::vector<int> Ls = {3, 4, 5, 6};
  for (auto L : Ls)
    {
      printf("TJModel: ALPS full spectrum test, chain N=%d\n", L);
      std::tie(bondlist, couplings) = tJchain(L, 1.0, 1.0);
      std::stringstream ss;
      ss << "data/tjfullspectrum/spectrum.chain." << L
	 << ".txt";	
      test_tjmodel_fullspectrum(bondlist, couplings, true,  ss.str());
    }  // Chains of length 3,4,5,6

  // Square 2x2
  printf("TJModel: ALPS full spectrum test, square 2x2\n");
  std::tie(bondlist, couplings) = square2x2(1.0, 1.0);
  test_tjmodel_fullspectrum(bondlist, couplings, true,
		    "data/tjfullspectrum/spectrum.square.2.txt");

  // Square 3x3
  printf("TJModel: ALPS full spectrum test, square 3x3\n");
  std::tie(bondlist, couplings) = square3x3(1.0, 1.0);
  test_tjmodel_fullspectrum(bondlist, couplings, true,
		    "data/tjfullspectrum/spectrum.square.3.txt");
   
  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of random all-to-all interactions 
  // by comparing to MATLAB results
  printf("TJModel: MATLAB full spectrum test, chain N=3\n");
  std::tie(bondlist, couplings) = randomAlltoAll3();
  test_tjmodel_fullspectrum(bondlist, couplings, false,
		    "data/tjfullspectrum/spectrum.allToAll.3.txt");


  printf("TJModel: MATLAB full spectrum test, chain N=4\n");
  std::tie(bondlist, couplings) = randomAlltoAll4();
  test_tjmodel_fullspectrum(bondlist, couplings, false,
		    "data/tjfullspectrum/spectrum.allToAll.4.txt");

  // test if complex tJ chain gives Hermitian matrix
  Ls = {3, 4, 5, 6};
  for (auto L : Ls)
    {
      printf("TJModel: complex hermitecity test, chain N=%d\n", L);
      std::tie(bondlist, couplings) = tJchain(L, 1.0, 1.0);
      couplings["T"] = complex(1.0, 1.0);
      test_tjmodel_hermitian(bondlist, couplings);
    }

  // test if complex tJ all-to-all gives Hermitian matrix
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("TJModel: complex hermitecity test, all-to-all N=%d\n", n_sites);
      std::tie(bondlist, couplings) = freefermion_alltoall_complex(n_sites);
      test_tjmodel_hermitian(bondlist, couplings);
    }
}

