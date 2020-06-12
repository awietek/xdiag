#include "catch.hpp"

#include <fstream>
#include <lila/all.h>
#include <hydra/all.h>
#include <random>

#include "testcases_hubbardmodel.h"

using namespace hydra::all;
using namespace lila;


// Test by comparing to exactly known solution for given quantumnumber
template <class coeff_t=double>
void test_hubbardmodel_e0(hydra::operators::BondList bondlist, 
		     hydra::operators::Couplings couplings, 
		     hydra::hilbertspaces::hubbard_qn qn, 
		     double e0)
{
  auto model = HubbardModel<coeff_t>(bondlist, couplings, qn);
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  auto eigs = lila::EigenvaluesSym(H);
  // printf("tje0 t: %f, J: %f, nup: %d, ndn: %d, eigs(0): %f, e0 %f\n",
  // 	 lila::real(couplings["T"]), lila::real(couplings["J"]),
  // 	 qn.n_upspins, qn.n_downspins, eigs(0), e0);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}


// Test by comparing to full spectrum of alps
void test_hubbardmodel_fullspectrum(hydra::operators::BondList bondlist, 
		       hydra::operators::Couplings couplings,
		       std::string filename)
{
  int n_sites = bondlist.n_sites();

  // Compute full spectrum in hydra
  lila::Vector<double> all_eigs;
  for (int nup=0; nup<=n_sites; ++nup)
    for (int ndn=0; ndn<=n_sites; ++ndn)
      {
	hubbard_qn qn = {nup, ndn};
	auto model = HubbardModel<double>(bondlist, couplings, qn);

	// Run Full ED
	auto H = model.matrix();

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
  for (int i=0; i<all_eigs.size(); ++i){
    REQUIRE(close(all_eigs(i), alps_eigs(i)));
  }
}


// Test if complex matrix is hermitian
void test_hubbardmodel_hermitian(hydra::operators::BondList bondlist, 
			    hydra::operators::Couplings couplings)
{
  int n_sites = bondlist.n_sites();
  for (int nup=0; nup<=n_sites; ++nup)
    for (int ndn=0; ndn<=n_sites-nup; ++ndn)
      {
 	auto model = HubbardModel<complex>(bondlist, couplings, {nup, ndn});
	auto H = model.matrix();
	REQUIRE(lila::close(H, lila::Herm(H)));
      }
}

TEST_CASE( "HubbardModel", "[HubbardModel]" )
{
  using namespace hydra::hubbardtestcases;
  BondList bondlist;
  Couplings couplings;

/*  
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("HubbardModel: free fermion random all-to-all test, N=%d\n", n_sites);
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
	    test_hubbardmodel_e0(bondlist, couplings, {nup,0}, e0);
	    test_hubbardmodel_e0(bondlist, couplings, {0,nup}, e0);
	  }
    }
  

  ///////////////////////////////////////////////////
  // Test Fermion all to all, free fermions (complex)
  for (int n_sites = 3; n_sites < 8; ++n_sites)
    {
      printf("HubbardModel: free fermion random all-to-all test (cplx), N=%d\n", n_sites);
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
	    test_hubbardmodel_e0<complex>(bondlist, couplings, {nup,0}, e0);
	    test_hubbardmodel_e0<complex>(bondlist, couplings, {0,nup}, e0);
	  }
    }

 */  
  /////////////////////////////////////////////////////////////////
  // Test of full spectrum of random all-to-all interactions 
  // by comparing to MATLAB results
  /*
  printf("TJModel: MATLAB full spectrum test, chain N=3\n");
  std::tie(bondlist, couplings) = randomAlltoAll3();
  test_hubbardmodel_fullspectrum(bondlist, couplings,
		    "data/hubbardfullspectrum/spectrum.allToAll.3.txt");
*/

  printf("HubbardModel: MATLAB full spectrum test, chain N=4, no Hubbard repulsion\n");
  std::tie(bondlist, couplings) = randomAlltoAll4NoU();
  test_hubbardmodel_fullspectrum(bondlist, couplings, 
		    "data/hubbardfullspectrum/spectrum.allToAll.N.4.U.0.txt");

  printf("HubbardModel: MATLAB full spectrum test, chain N=4, Hubbard repulsion\n");
  std::tie(bondlist, couplings) = randomAlltoAll4();
  test_hubbardmodel_fullspectrum(bondlist, couplings, 
		    "data/hubbardfullspectrum/spectrum.allToAll.N.4.U.5.txt");

}

