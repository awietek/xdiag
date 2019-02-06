#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <limits>

#include <lila/all.h>
#include <hydra/all.h>

#include "spinlessfermioned.options.h"

int main(int argc, char* argv[])
{
  using hydra::models::SpinlessFermions;
  using hydra::utils::range;
  using hydra::operators::BondList;
  using hydra::operators::read_bondlist;
  using hydra::symmetries::read_charactertable;

  double t = 1.;
  double V = 1.;
  std::string outfile;
  std::string latticefile;
  int np = -1;
  std::string representation;
  parse_cmdline(t, V, outfile, latticefile, np, representation, argc, argv);

  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  std::cerr << "HeisenbergThermo Error in opening outfile: " 
		    << "Could not open file with filename ["
		    << outfile << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }

  // Create model from file
  BondList bondlist = read_bondlist(latticefile);
  int n_sites = bondlist.n_sites();
  printf("n_sites: %d\n", n_sites);
  BondList hopping_list = bondlist.bonds_of_type("HOP");
  std::vector<std::pair<int, int>> hoppings;
  for (auto bond : hopping_list)
    hoppings.push_back({bond.sites()[0], bond.sites()[1]});

  BondList interaction_list = bondlist.bonds_of_type("INT");
  std::vector<std::pair<int, int>> interactions;
  for (auto bond : hopping_list)
    interactions.push_back({bond.sites()[0], bond.sites()[1]});

  auto model = SpinlessFermions(n_sites, hoppings, interactions);

  // Select qns to compute
  std::vector<int> qns;
  if (np == -1) 
    {
      for (auto qn : range<int>(0, n_sites + 1))
	qns.push_back(qn);
    }
  else
    qns = {np};

  // Compute energies
  int qn_idx = 0;
  for (auto qn : qns)
    {
      if (representation == "")
	{
	  printf("No lattice symmetry used\n");
	  printf("Creating spinless fermion Hamiltonian for n_particles=%d...\n", qn);
	  auto hamilton = model.matrix(t, V, qn);
	  printf("dim: %d\n", hamilton.nrows());
	  printf("Done\n");
	  
	  if (hamilton.nrows()>0)
	    {
	      printf("Computing eigenvalues ...\n");
	      auto eigs = lila::EigenvaluesH(hamilton);
	      double e0 = eigs(0); 
	      printf("e0: %.20g\n", e0);
	      printf("Done\n");
	    }

	}  
      else
	{
	  printf("Using symmetry sector %s\n", representation.c_str());
	  auto character_table = read_charactertable(latticefile);
	  printf("Creating spinless fermion Hamiltonian for n_particles=%d...\n", qn);
	  auto hamilton = model.matrix(t, V, qn, character_table, representation);
	  // LilaPrint(hamilton);
	  printf("dim: %d\n", hamilton.nrows());
	  printf("Done\n");
	  assert(lila::close(Conj(Transpose(hamilton)), hamilton));

	  if (hamilton.nrows()>0)
	    {
	      printf("Computing eigenvalues ...\n");
	      auto eigs = lila::EigenvaluesH(hamilton);
	      double e0 = eigs(0); 
	      printf("e0: %.20g\n", e0);
	      printf("Done\n");
	    }
	}  
      ++qn_idx;
    }

  // // Write to outfile
  // std::stringstream line;
  // double total_e0 = *std::min_element(e0s.begin(), e0s.end());
  // line << "# gs energy (ED): " << total_e0 << "\n";
  
  return EXIT_SUCCESS;
}
