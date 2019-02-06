#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <limits>

#include <lila/all.h>
#include <hydra/all.h>

#include "heisenberged.options.h"

int main(int argc, char* argv[])
{
  using hydra::models::HeisenbergModel;
  using hydra::utils::range;
  using hydra::operators::BondList;
  using hydra::operators::read_bondlist;
  using hydra::symmetries::read_charactertable;

  double J = 1.;
  std::string outfile;
  std::string latticefile;
  int nup = -1;
  std::string representation;
  parse_cmdline(J, outfile, latticefile, nup, representation, argc, argv);

  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  std::cerr << "HeisenbergED Error in opening outfile: " 
		    << "Could not open file with filename ["
		    << outfile << "] given. Abort." << std::endl;
	  exit(EXIT_FAILURE);
	}
    }

  // Create model from file
  BondList bondlist = read_bondlist(latticefile);
  int n_sites = bondlist.n_sites();
  printf("n_sites: %d\n", n_sites);
  BondList hopping_list = bondlist.bonds_of_type("HB");
  std::vector<std::pair<int, int>> hoppings;
  for (auto bond : hopping_list)
    hoppings.push_back({bond.sites()[0], bond.sites()[1]});

  auto model = HeisenbergModel(n_sites, hoppings);

  // Select qns to compute
  std::vector<int> qns;
  if (nup == -1) 
    {
      for (auto qn : range<int>(0, n_sites + 1))
	qns.push_back(qn);
    }
  else
    qns = {nup};

  // Compute energies
  int qn_idx = 0;
  for (auto qn : qns)
    {
      if (representation == "")
	{
	  printf("No lattice symmetry used\n");
	  printf("Creating Heisenberg Hamiltonian for n_particles=%d...\n", qn);
	  auto hamilton = model.matrix(J, qn);
	  printf("dim: %d\n", hamilton.nrows());
	  printf("Done\n");
	  

	  if (hamilton.nrows()>0)
	    {
	      printf("Computing eigenvalues ...\n");
	      auto eigs = lila::EigenvaluesH(hamilton);
	      double e0 = eigs(0); 
	      printf("e0: %f\n", e0);
	      printf("Done\n");
	    }

	}  
      else
	{
	  printf("Using symmetry sector %s\n", representation.c_str());
	  auto character_table = read_charactertable(latticefile);
	  printf("Creating Heisenberg Hamiltonian for n_particles=%d...\n", qn);
	  auto hamilton = model.matrix(J, qn, character_table, representation);
	  // LilaPrint(hamilton);
	  printf("dim: %d\n", hamilton.nrows());
	  printf("Done\n");
	  assert(lila::close(Conj(Transpose(hamilton)), hamilton));

	  if (hamilton.nrows()>0)
	    {
	      printf("Computing eigenvalues ...\n");
	      auto eigs = lila::EigenvaluesH(hamilton);
	      double e0 = eigs(0); 
	      printf("e0: %f\n", e0);
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
