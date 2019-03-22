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
  using namespace hydra::utils;
  using namespace hydra::operators;
  using namespace hydra::thermodynamics;

  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string evaluate = "evaluate";
  std::string method = "exact";
  std::string temperaturefile;
  int seed = 42;
  double precision = 1e-12;
  int neval = 1;
  int iters = 1000;
  int nup = -1;
  int ndown = -1;
  int np = -1;
  bool loseevals = false;

  parse_cmdline(outfile, latticefile, couplingfile, evaluate, method, temperaturefile, 
		seed, precision, neval, iters, nup, ndown, np, loseevals, argc, argv);

  // Input checks
  check_if_file_exists(latticefile);
  check_if_file_exists(couplingfile);
  check_if_contained_in(method, {"exact", "TPQ"});
  check_if_contained_in(evaluate, {"evaluate", "notevaluate", "onlyevaluate"});
  if (!(evaluate == "notevaluate"))
    check_if_file_exists(temperaturefile);
  assert(neval >= 0);

  // open outfile
  std::ofstream of;
  of.open(outfile, std::ios::out);
  check_if_file_exists(outfile);
  of << "### method: " << method << "\n";

  // Parse bondlist and couplings from file
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
    printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), std::imag(c.second));
  int n_sites = bondlist.n_sites();

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
  if (evaluate != "notevaluate")
    temperatures = lila::ReadVector<double>(temperaturefile);
  

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
      auto model = HubbardModel(bondlist, couplings, qn); 
      printf("dim: %s\n", FormatWithCommas(model.dim()).c_str());

      if (method=="exact")
	{
	  auto hamilton = model.matrix();
	  printf("Done\n");

	  printf("Computing eigenvalues ...\n");
	  auto eigs = lila::EigenvaluesH(hamilton);
	  double e0 = eigs(0); 
	  printf("e0: %f\n", e0);
	  e0s.push_back(e0);
	  printf("Done\n");
      
	  // Write eigenvalues to file
	  if (!loseevals)
	    {
	      of << "## BLOCK: nup=" << qn.n_upspins << ", ndown=" << qn.n_downspins
		 << "\n";
	      of << "# eigenvalues\n";
	      for (auto e : eigs)
		of << std::setprecision(20) <<  e << "\n";
	    }

	  // Compute thermodynamic quantities for all temperatures
	  if (evaluate != "notevaluate")
	    {
	      if (loseevals)
		of << "## BLOCK: nup=" << qn.n_upspins << ", ndown=" << qn.n_downspins
		   << "\n";

	      of << "# temperature partition energy quadmoment\n";
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
		  auto partition = lila::Sum(exp_eigs);
		  auto energy = lila::Dot(eigs, exp_eigs);
		  auto quadmoment = lila::Dot(sq_eigs, exp_eigs);
		  partitions_for_qn[t_idx][qn_idx] = partition;
		  energies_for_qn[t_idx][qn_idx] = energy;
		  quad_moments_for_qn[t_idx][qn_idx] = quadmoment;
		  of << std::setprecision(20)  
		     << t << " " << partition << " " 
		     << energy / partition << " " 
		     << quadmoment / partition << "\n";
		  ++t_idx;
		}
	    }
	}
      else if (method=="TPQ")
	{
	  printf("Done\n");
	  auto multiply = [&model](const Vector<double>& v, Vector<double>& w) {
	    static int iter=0;
	    auto t1 = Clock::now();
	    model.apply_hamiltonian(v, w);
	    auto t2 = Clock::now();
	    printf("iter: %d, time MVM: %3.4f\n", iter, secs(t2-t1).count()); 
	    ++iter;
	  };
  
	  int64 dim = model.dim();
	  auto lzs = Lanczos<double, decltype(multiply)>
	    (dim, seed + 1234321*qn_idx, iters, precision, neval, multiply);
	  Vector<double> eigs = lzs.eigenvalues();
	  printf("lzs e %20.18g\n", eigs(0));
	  printf("Done\n");

	  // Write alphas and betas to file
	  if (!loseevals)
	    {
	      of << "## BLOCK: nup=" << qn.n_upspins << ", ndown=" << qn.n_downspins
		 << "\n";
	      of << "# alphas betas\n";
	      for (int i=0; i<lzs.tmatrix().diag().size(); ++i)
		of << std::setprecision(20) 
		   << lzs.tmatrix().diag()(i) << " " 
		   << lzs.tmatrix().offdiag()(i)<< "\n";
	    }

	  // Compute thermodynamic quantities for all temperatures
	  if (evaluate != "notevaluate")
	    {
	      if (loseevals)
		of << "## BLOCK: nup=" << qn.n_upspins << ", ndown=" << qn.n_downspins
		   << "\n";

	      of << "# temperature partition energy quadmoment\n";
	      int t_idx = 0;
	      auto tmat = lzs.tmatrix();
	      auto teig = Eigen(tmat);
	      auto Q = teig.eigenvectors;
	      auto teigs = teig.eigenvalues;
	      double shift = *std::min_element(eigs.begin(), eigs.end());
	      for (double t : temperatures)
		{
		  double beta = 1. / t;
		  auto diag = Zeros(tmat.size(), tmat.size());
		  for (int j = 0; j < tmat.size(); ++j)
		    diag(j, j) = exp(-(beta / 2.) * (teigs(j) - shift));
		  auto B = Mult(Q, diag)
		  
		  auto exp_eigs = eigs;
		  lila::Map(exp_eigs, 
			    [beta, e0](double& e) { 
			      e = std::exp(-beta * (e - e0)); 
			    });
		  auto sq_eigs = eigs;
		  lila::Map(sq_eigs, [](double& e) { e = e*e; });
		  auto partition = lila::Sum(exp_eigs);
		  auto energy = lila::Dot(eigs, exp_eigs);
		  auto quadmoment = lila::Dot(sq_eigs, exp_eigs);
		  partitions_for_qn[t_idx][qn_idx] = partition;
		  energies_for_qn[t_idx][qn_idx] = energy;
		  quad_moments_for_qn[t_idx][qn_idx] = quadmoment;
		  of << std::setprecision(20)  
		     << t << " " << partition << " " 
		     << energy / partition << " " 
		     << quadmoment / partition << "\n";
		  ++t_idx;
		}
	    }

	  
	}
      ++qn_idx;
    }

  // Combine quantities from different quantum numbers
  if (evaluate != "notevaluate")
    { 
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
      of << "## COMBINATION\n";
      std::stringstream line;
      double total_e0 = *std::min_element(e0s.begin(), e0s.end());
      of << "# gs energy (ED): " << total_e0 << "\n";
      of << "# temperature    partition   energy    quadmoments      specificheats\n";
      t_idx = 0;
      for (double t : temperatures)
	{
	  of << t << " " << partitions[t_idx] << " " << energies[t_idx] << " " 
	     << quad_moments[t_idx] << " " << specific_heats[t_idx] <<  "\n";
	  ++t_idx;
	}
    }
  return EXIT_SUCCESS;
}
