#include "catch.hpp"

#include <fstream>
#include <lila/all.h>
#include <hydra/all.h>

using namespace hydra::all;
using namespace lila;

void test_tjmodel(hydra::operators::BondList bondlist, 
		  hydra::operators::Couplings couplings, 
		  hydra::hilbertspaces::hubbard_qn qn, 
		  double e0)
{
  auto model = TJModel<double>(bondlist, couplings, qn);

  // Run Lanczos
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  auto eigs = lila::EigenvaluesSym(H);
  REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}



void test_tjmodel_alps(hydra::operators::BondList bondlist, 
		       hydra::operators::Couplings couplings,
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
	bool ninj_term = true;
	auto H = model.matrix(ninj_term);
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

  REQUIRE(all_eigs.size() == alps_eigs.size());
  for (int i=0; i<all_eigs.size(); ++i)
    REQUIRE(close(all_eigs(i), alps_eigs(i)));
  
}


void test_tjmodel_complex(hydra::operators::BondList bondlist, 
			  hydra::operators::Couplings couplings, 
			  hydra::hilbertspaces::hubbard_qn qn)
{
  auto model = TJModel<complex>(bondlist, couplings, qn);
  auto H = model.matrix();
  REQUIRE(lila::close(H, lila::Herm(H)));
  
  // auto eigs = lila::EigenvaluesSym(H);
  // REQUIRE(std::abs(e0 - eigs(0)) < 1e-6);
}

TEST_CASE( "TJModel test", "[TJModel]" ) {

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
    	test_tjmodel(bondlist, couplings, qns[i], e0s[i]);
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
	test_tjmodel(bondlist, couplings, qns[i], e0s[i]);
    }

  }


  // compare full spectrum of chain with alps
  {
    double t=1.0;
    double J=1.0;
    Couplings couplings;
    couplings["T"] = t;
    couplings["J"] = J;

    // Chains of length 3,4,5,6
    {
      std::vector<int> Ls = {3, 4, 5, 6};
      for (auto L : Ls)
	{
	  BondList bondlist;
	  for (int s=0; s<L; ++s)
	    {
	      bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % L});
	      bondlist << Bond("HEISENBERG", "J", {s, (s+1) % L});
	    }
	  std::stringstream ss;
	  ss << "data/tjfullspectrum/spectrum.chain." << L
	     << ".txt";	
	  test_tjmodel_alps(bondlist, couplings, ss.str());
	}
    }  // Chains of length 3,4,5,6

    // Square 2x2
    {
      BondList bondlist;
      bondlist << Bond("HUBBARDHOP", "T", {0, 1});
      bondlist << Bond("HUBBARDHOP", "T", {1, 0});
      bondlist << Bond("HUBBARDHOP", "T", {2, 3});
      bondlist << Bond("HUBBARDHOP", "T", {3, 2});
      bondlist << Bond("HUBBARDHOP", "T", {0, 2});
      bondlist << Bond("HUBBARDHOP", "T", {2, 0});
      bondlist << Bond("HUBBARDHOP", "T", {1, 3});
      bondlist << Bond("HUBBARDHOP", "T", {3, 1});
      bondlist << Bond("HEISENBERG", "J", {0, 1});
      bondlist << Bond("HEISENBERG", "J", {1, 0});
      bondlist << Bond("HEISENBERG", "J", {2, 3});
      bondlist << Bond("HEISENBERG", "J", {3, 2});
      bondlist << Bond("HEISENBERG", "J", {0, 2});
      bondlist << Bond("HEISENBERG", "J", {2, 0});
      bondlist << Bond("HEISENBERG", "J", {1, 3});
      bondlist << Bond("HEISENBERG", "J", {3, 1});
      test_tjmodel_alps(bondlist, couplings,
			"data/tjfullspectrum/spectrum.square.2.txt");
    }  // Square 2x2

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
      
      test_tjmodel_alps(bondlist, couplings,
			"data/tjfullspectrum/spectrum.square.3.txt");
    }  // Square 3x3

  }  // compare full spectrum of chain with alps


  // test if complex tJ model gives Hermitian matrix
  {  
    Couplings couplings;
    couplings["T"] = complex(1.0, 1.0);
    couplings["J"] = 1.0;

    std::vector<int> Ls = {3, 4, 5, 6};
    for (auto L : Ls)
      {
	BondList bondlist;
	for (int s=0; s<L; ++s)
	  {
	    bondlist << Bond("HUBBARDHOP", "T", {s, (s+1) % L});
	    bondlist << Bond("HEISENBERG", "J", {s, (s+1) % L});
	  }

	std::vector<hubbard_qn> qns; 
	for (int nup=0; nup<=L; ++nup)
	  for (int ndn=0; ndn<L - nup; ++ndn)
	    {
	      hubbard_qn qn = {nup, ndn};
	      test_tjmodel_complex(bondlist, couplings, qn);
	    }
      }

  }

}
