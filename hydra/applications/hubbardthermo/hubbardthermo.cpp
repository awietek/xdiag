#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>

#include <lila/all.h>
#include <hydra/all.h>

#include "hubbardthermo.options.h"

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::models::HubbardModel;
  using hydra::utils::range;
  using hydra::operators::BondList;
  using hydra::operators::read_bondlist;
  using hydra::thermodynamics::combined_quantity;

  double T = 0.1;
  double t = 1.;
  double U = 4.;
  std::string outfile;
  std::string latticefile;
  std::string temperaturefile;
  int nup = -1;
  int ndown = -1;
  int np = -1;

  parse_cmdline(T, t, U, outfile, latticefile, temperaturefile, nup, ndown, 
		np, argc, argv);

  // Open outfile
  std::ofstream of;
  of.open(outfile, std::ios::out);
  if(of.fail()) 
    {
      std::cerr << "HubbardThermo Error in opening outfile: " 
		<< "Could not open file with filename ["
		<< outfile << "] given. Abort." << std::endl;
      exit(EXIT_FAILURE);
    }


  // Parse hoppings from file
  BondList bondlist = read_bondlist(latticefile);
  int n_sites = bondlist.n_sites();
  printf("n_sites: %d\n", n_sites);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  std::vector<std::pair<int, int>> hoppings;
  for (auto bond : hopping_list)
    hoppings.push_back({bond.sites()[0], bond.sites()[1]});
  auto model = HubbardModel(n_sites, hoppings);


  // Select qns to compute
  std::vector<hubbard_qn> qns;
  if ((nup == -1)  || (ndown == -1))
    {
      for (auto up : range<int>(0, n_sites + 1))
	for (auto down : range<int>(0, n_sites + 1))
	  if ((np==-1) || (up + down == np))
	    qns.push_back({up, down});
    }
  else qns = {{nup,ndown}};


  // Set temperatures to compute
  lila::Vector<double> temperatures;
  if (temperaturefile == "") 
    {
      temperatures.resize(1);
      temperatures(0) = T;
    }
  else temperatures = lila::ReadVector<double>(temperaturefile);

  // Allocate result data
  std::vector<std::vector<double>> partitions_for_qn(temperatures.size());
  std::vector<std::vector<double>> energies_for_qn(temperatures.size());
  std::vector<std::vector<double>> quad_moments_for_qn(temperatures.size());
  for (int i=0; i < (int)temperatures.size(); ++i)
    {
      partitions_for_qn[i].resize(qns.size());
      energies_for_qn[i].resize(qns.size());
      quad_moments_for_qn[i].resize(qns.size());
    }

  // Compute thermodymanic quantities for all quantum numbers
  std::vector<double> e0s;
  int qn_idx = 0;
  for (auto qn : qns)
    {
      printf("Creating Hubbard Hamiltonian for n_upspins=%d, n_downspins=%d...\n", 
	     qn.n_upspins, qn.n_downspins);
      auto hamilton = model.matrix(t, U, qn);
      printf("Done\n");

      printf("Computing eigenvalues ...\n");
      auto eigs = lila::EigenvaluesH(hamilton);
      double e0 = eigs(0); 
      printf("e %f %f %f %f\n", e0, eigs(1), t, U);

      e0s.push_back(e0);
      printf("Done\n");
  
      // Compute thermodynamic quantities for all temperatures
      int t_idx = 0;
      for (double t : temperatures)
	{
	  double beta = 1. / t;
	  auto exp_eigs = eigs;
	  lila::Map(exp_eigs, 
		    [beta, e0](double& e) { 
		      e = std::exp(-beta * (e - e0)); 
		    });
	  auto sq_eigs = eigs;
	  lila::Map(sq_eigs, [](double& e) { e = e*e; });
	  partitions_for_qn[t_idx][qn_idx] = lila::Sum(exp_eigs);
	  energies_for_qn[t_idx][qn_idx] = lila::Dot(eigs, exp_eigs);
	  quad_moments_for_qn[t_idx][qn_idx] = lila::Dot(sq_eigs, exp_eigs);
	  printf("%f %f\n", e0, energies_for_qn[t_idx][qn_idx] / partitions_for_qn[t_idx][qn_idx]);

	  ++t_idx;
	}
      ++qn_idx;
    }

  // Combine different quantum numbers and dump
  std::vector<double> partitions(temperatures.size(), 0.);
  std::vector<double> energies(temperatures.size(), 0.);
  std::vector<double> quad_moments(temperatures.size(), 0.);
  std::vector<double> specific_heats(temperatures.size(), 0.);

  int t_idx = 0;
  for (double t : temperatures)
    {
      double beta = 1. / t;
      partitions[t_idx] = combined_quantity(partitions_for_qn[t_idx], e0s, beta);
      energies[t_idx] = combined_quantity(energies_for_qn[t_idx], e0s, beta);
      quad_moments[t_idx] = combined_quantity(quad_moments_for_qn[t_idx], e0s, beta);

      energies[t_idx] /= partitions[t_idx];
      quad_moments[t_idx] /= partitions[t_idx];
      specific_heats[t_idx] = beta * beta * 
	(quad_moments[t_idx] - energies[t_idx] * energies[t_idx]);
      ++t_idx;
    }

  // Write to outfile
  std::stringstream line;
  double total_e0 = *std::min_element(e0s.begin(), e0s.end());
  line << "# gs energy (ED): " << total_e0 << "\n";
  of << line.str();
  line.str("");
  line << "# temperature    partition   energy    specific heat\n";
  std::cout << line.str();
  of << line.str();
  line.str("");

  t_idx = 0;
  for (double t : temperatures)
    {
      double beta = 1. / t;
      line << t << " " << partitions[t_idx] << " " << energies[t_idx] << " " 
	   << specific_heats[t_idx] << "\n";
      std::cout << line.str();
      of << line.str();
      line.str("");
      ++t_idx;
    }
  
  return EXIT_SUCCESS;
}
