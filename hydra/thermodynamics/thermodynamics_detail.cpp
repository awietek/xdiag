#include <algorithm>
#include <cmath>
#include <hydra/thermodynamics/thermodynamics_detail.h>

namespace hydra { namespace thermodynamics { namespace detail {

      void combine_thermodynamics
      (const std::vector<double> temperatures,
       const std::vector<double>& e0s,
       std::vector<std::vector<double>>& partitions_for_qn,
       std::vector<std::vector<double>>& energies_for_qn,
       std::vector<std::vector<double>>& quad_moments_for_qn,
       std::vector<double>& partitions,
       std::vector<double>& energies,
       std::vector<double>& quad_moments,
       std::vector<double>& specific_heats)
      {
	int n_temperatures = temperatures.size();
	assert((int)partitions_for_qn.size() == n_temperatures);
	assert((int)energies_for_qn.size() == n_temperatures);
	assert((int)quad_moments_for_qn.size() == n_temperatures);

	partitions.resize(n_temperatures, 0);
	energies.resize(n_temperatures, 0);
	quad_moments.resize(n_temperatures, 0);
	specific_heats.resize(n_temperatures, 0);

	int t_idx = 0;
	for (double t : temperatures)
	  {
	    double beta = 1. / t;
	    partitions[t_idx] = combined_quantity(partitions_for_qn[t_idx], 
						  e0s, beta);
	    energies[t_idx] = combined_quantity(energies_for_qn[t_idx], 
						e0s, beta);
	    quad_moments[t_idx] = combined_quantity(quad_moments_for_qn[t_idx], 
						    e0s, beta);

	    energies[t_idx] /= partitions[t_idx];
	    quad_moments[t_idx] /= partitions[t_idx];

	    specific_heats[t_idx] = beta * beta * 
	      (quad_moments[t_idx] - energies[t_idx] * energies[t_idx]);

	    for (int qn_idx=0; qn_idx < (int)partitions_for_qn[t_idx].size(); ++qn_idx)
	      {
		energies_for_qn[t_idx][qn_idx] /= 
		  partitions_for_qn[t_idx][qn_idx];
		quad_moments_for_qn[t_idx][qn_idx] /= 
		  partitions_for_qn[t_idx][qn_idx];	      
	      }
	    ++t_idx;
	  }
      }

      double combined_quantity(const std::vector<double>& quantities, 
			       const std::vector<double>& e0s, double beta)
      {
	assert(quantities.size() == e0s.size());
	double total_e0 = *std::min_element(e0s.begin(), e0s.end());
	double combined_quantity = 0.;
	for (int idx = 0; idx < (int)e0s.size(); ++idx)
	  {
	    double pre = std::exp(-beta * (e0s[idx] - total_e0));
	    combined_quantity += pre * quantities[idx];
	  }
	return combined_quantity;
      }


    }
  }
}
