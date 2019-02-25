#include "continuedfraction.h"

namespace hydra { namespace dynamics {

    std::complex<double> continued_fraction(const std::complex<double> z, 
					    const std::vector<double>& alphas,
					    const std::vector<double>& betas, 
					    int level)
    {
      if (level==(int)alphas.size())
	return 0.0;

      else if(level == 0)
	return 1./ (z - alphas[0] - continued_fraction(z, alphas, betas, 1));
      else
        return (betas[level-1]*betas[level-1]) /			\
            (z - alphas[level] - continued_fraction(z, alphas, betas, level+1));
    }
    

  }
}
