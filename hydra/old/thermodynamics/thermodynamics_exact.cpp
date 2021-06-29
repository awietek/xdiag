#include <cassert>
#include <cmath>

#include <hydra/thermodynamics/thermodynamics_exact.h>
#include <hydra/thermodynamics/thermodynamics_detail.h>
#include <hydra/utils/iochecks.h>
#include <hydra/models/hubbardmodel.h>

namespace hydra { namespace thermodynamics {

    template <class model_t>
    thermodynamics_exact_result_t
    thermodynamics_exact(const BondList& bondlist, 
			 const Couplings& couplings,
			 const std::vector<typename model_t::qn_t>& qns,
			 const std::vector<double>& temperatures,
			 bool keep_evs)
    {
      // Allocate result data
      int n_temperatures = temperatures.size();
      int n_qns = qns.size();
      thermodynamics_exact_result_t result;
      result.partitions_for_qn.resize(n_temperatures);
      result.energies_for_qn.resize(n_temperatures);
      result.quad_moments_for_qn.resize(n_temperatures);
      for (int i=0; i < (int)n_temperatures; ++i)
	{
	  result.partitions_for_qn[i].resize(n_qns, 0);
	  result.energies_for_qn[i].resize(n_qns, 0);
	  result.quad_moments_for_qn[i].resize(n_qns, 0);
	}
      result.eigenvalues.resize(n_qns);

      // Loop over all quantum numbers for individual thermodynamics
      int qn_idx = 0;
      for (auto qn : qns)
	{
	  
	  // Create and diagonalize Hamiltonian
	  auto model = model_t(bondlist, couplings, qn); 
	  auto hamilton = model.matrix();
	  auto eigs = lila::EigenvaluesSym(hamilton);
	  double e0 = eigs(0);
	  result.e0s.push_back(e0);

	  // Evaluate energy and quad moments at every temperature
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
	      double partition = lila::Sum(exp_eigs);
	      double energy = lila::Dot(eigs, exp_eigs);
	      double quadmoment = lila::Dot(sq_eigs, exp_eigs);
	      result.partitions_for_qn[t_idx][qn_idx] = partition;
	      result.energies_for_qn[t_idx][qn_idx] = energy;
	      result.quad_moments_for_qn[t_idx][qn_idx] = quadmoment;
	      ++t_idx;
	    }
	  if (keep_evs) result.eigenvalues[qn_idx] = eigs;
	  ++qn_idx;
	}  // for (auto qn : qns)

      // Combine results from different quantum numbers
      detail::combine_thermodynamics
	(temperatures, result.e0s, result.partitions_for_qn,
	 result.energies_for_qn, result.quad_moments_for_qn,
	 result.partitions, result.energies, 
	 result.quad_moments, result.specific_heats);
      return result;
    }

    
    template <class qn_t>
    void write_thermodynamics_exact
    (const std::vector<qn_t>& qns, 
     const std::vector<double>& temperatures,
     const thermodynamics_exact_result_t& result,
     std::string outfile, bool write_evs)
    {
      std::ofstream of;
      of.open(outfile, std::ios::out);
      utils::check_if_file_exists(outfile);
      of << "### method: exact\n";
      int qn_idx = 0;
      for (auto qn : qns)
	{
	  of << std::setprecision(20) << "## BLOCK: " << qn << "\n"	      
	     << "# gs energy: " << result.e0s[qn_idx] << "\n";
	  if (write_evs)
	    {
	      of << "# eigenvalues\n";
	      for (auto e : result.eigenvalues[qn_idx])
		of << std::setprecision(20) <<  e << "\n";
	    }
	  of << "# temperature partition energy quadmoment\n";
	  int t_idx = 0;
	  for (double t : temperatures)
	    {
	      of << std::setprecision(20)  
		 << t << " " << result.partitions_for_qn[t_idx][qn_idx] << " " 
		 << result.energies_for_qn[t_idx][qn_idx] << " " 
		 << result.quad_moments_for_qn[t_idx][qn_idx] << "\n";
	      ++t_idx;
	    }
	  ++qn_idx;
	}
      of << "## COMBINATION\n";
      double total_e0 = *std::min_element(result.e0s.begin(), result.e0s.end());
      of << "# gs energy: " << std::setprecision(20) << total_e0 << "\n";
      of << "# temperature  partition  energy  quadmoments  specificheats\n";

      int t_idx = 0;
      for (double t : temperatures)
	{
	  of << std::setprecision(20)  
	     << t << " " 
	     << result.partitions[t_idx] << " " 
	     << result.energies[t_idx] << " " 
	     << result.quad_moments[t_idx] << " " 
	     << result.specific_heats[t_idx] << "\n";
	  ++t_idx;
	}
    }
    
    // Hubbard model 
    using models::HubbardModel;  
    template thermodynamics_exact_result_t 
    thermodynamics_exact<HubbardModel<double>>
    (const BondList& bondlist, const Couplings& couplings,
     const std::vector<HubbardModel<double>::qn_t>& qns,
     const std::vector<double>& temperatures, bool keep_evs);


    template thermodynamics_exact_result_t 
    thermodynamics_exact<HubbardModel<complex>>
    (const BondList& bondlist, const Couplings& couplings,
     const std::vector<HubbardModel<complex>::qn_t>& qns,
     const std::vector<double>& temperatures, bool keep_evs);

    template 
    void write_thermodynamics_exact<HubbardModel<double>::qn_t>
    (const std::vector<HubbardModel<double>::qn_t>& qns, 
     const std::vector<double>& temperatures,
     const thermodynamics_exact_result_t& result,
     std::string outfile, bool write_evs);


  }
}
