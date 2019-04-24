#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <chrono>

#include <lila/all.h>
#include <hydra/all.h>

#include "hubbardthermo.options.h"

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::models::HubbardModel;
  using namespace hydra::utils;
  using namespace hydra::operators;
  using namespace hydra::thermodynamics;
  using namespace lila;
  using Clock = std::chrono::high_resolution_clock;
  using secs = std::chrono::duration<float>;

  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string method = "exact";
  std::string temperaturefile;
  int seed = 42;
  double precision = 1e-12;
  int neval = 1;
  int iters = 1000;
  int nup = -1;
  int ndown = -1;
  int np = -1;
  bool writeevals = false;
  parse_cmdline(outfile, latticefile, couplingfile, method, temperaturefile, 
		seed, precision, neval, iters, nup, ndown, np, writeevals, 
		argc, argv);

  // Input checks
  check_if_file_exists(latticefile);
  check_if_file_exists(couplingfile);
  check_if_contained_in(method, {"exact", "TPQ"});
  check_if_file_exists(temperaturefile);
  assert(neval >= 0);

  // Parse bondlist/couplings/temperatures from file
  BondList bondlist = read_bondlist(latticefile);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  for (auto bond : hopping_list)
    printf("hopping %s %d %d\n", bond.coupling().c_str(), 
	   bond.sites()[0], bond.sites()[1]);
  BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
  for (auto bond : interaction_list)
    printf("interaction %s %d %d\n", bond.coupling().c_str(), 
	   bond.sites()[0], bond.sites()[1]);
  if (couplingfile == "") couplingfile = latticefile;
  Couplings couplings = read_couplings(couplingfile);
  for (auto c : couplings)
    printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), 
	   std::imag(c.second));
  auto temperatures = lila::ReadVector<double>(temperaturefile);

  // Select qns to compute
  std::vector<hubbard_qn> qns;
  int n_sites = bondlist.n_sites();
  if ((nup == -1)  || (ndown == -1))
    {
      for (int up=0; up<n_sites + 1; ++up)
	for (int down=0; down<n_sites + 1; ++down)
	  if ((np==-1) || (up + down == np))
	    qns.push_back({up, down});
    }
  else qns = {{nup,ndown}};

  // Choose method to compute
  if (method=="exact")
    {
      auto result = thermodynamics_exact<HubbardModel>
	(bondlist, couplings, qns, temperatures);
      write_thermodynamics_exact(qns, temperatures, result, 
				 outfile, writeevals);
    }      
  else if (method=="TPQ")
    {
      auto result = thermodynamics_tpq<HubbardModel>
	(bondlist, couplings, qns, temperatures, seed, iters, precision, neval);
      write_thermodynamics_tpq(qns, temperatures, result, outfile);
    }

  return EXIT_SUCCESS;
}
