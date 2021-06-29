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

  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string ensemble = "canonical";
  std::string temperaturefile;
  int np = -1;
  double mu = 0.;
  parse_cmdline(outfile, latticefile, couplingfile, ensemble, temperaturefile, 
		np, mu, argc, argv);

  // Input checks
  check_if_file_exists(latticefile);
  check_if_file_exists(couplingfile);
  check_if_contained_in(ensemble, {"canonical", "grandcanonical"});
  check_if_file_exists(temperaturefile);

  // Parse bondlist/couplings/temperatures from file
  auto bondlist = read_bondlist(latticefile);
  auto couplings = read_couplings(couplingfile);
  auto temperatures = lila::ReadVector<double>(temperaturefile);

  // Select qns to compute
  int n_sites = bondlist.n_sites();
  std::vector<hubbard_qn> qns;
  if (ensemble == "canonical")
    {
      if (np < 0) np = n_sites;
      for (int nup=0; nup<=np; ++nup)
	{
	  int ndown = np - nup;
	  qns.push_back({nup, ndown});
	}
    }
  else if (ensemble == "grandcanonical")
    {
      for (int nup=0; nup<=n_sites; ++nup)
	for (int ndown=0; ndown<=n_sites; ++ndown)
	  qns.push_back({nup, ndown});
      
      couplings["MU"] = mu;
      for (int site=0; site<n_sites; ++site)
	bondlist << Bond("HUBBARDMU", "MU", {site});      
    }

  auto result = thermodynamics_exact<HubbardModel<double>>
    (bondlist, couplings, qns, temperatures);
  write_thermodynamics_exact(qns, temperatures, result, outfile);

  return EXIT_SUCCESS;
}
