#include <cassert>
#include <algorithm>

#include "thermodynamics.h"

namespace hydra { namespace thermodynamics {

    double combined_quantity(const std::vector<double>& quantities, 
			     const std::vector<double>& e0s, double beta)
    {
      assert(quantities.size() == e0s.size());
      double total_e0 = *std::min_element(e0s.begin(), e0s.end());
      double combined_quantity = 0.;
      for (int idx = 0; idx < e0s.size(); ++idx)
	{
	  double pre = std::exp(-beta * (e0s[idx] - total_e0));
	  combined_quantity += pre * quantities[idx];
	}
      return combined_quantity;
    }

  }
}
