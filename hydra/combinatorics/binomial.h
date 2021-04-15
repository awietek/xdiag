#ifndef HYDRA_COMBINATORICS_BINOMIAL_
#define HYDRA_COMBINATORICS_BINOMIAL_

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

namespace hydra {
namespace combinatorics {

  int64 binomial(int n, int k);

  constexpr int64 binom(int n, int k) {
    if (k > n || k < 0)
      return 0;
    int64 res = 1;
    for (int i = 1; i <= k; i++)
      res = (res * (n - i + 1)) / i;
    return res;
  } 

} // namespace combinatorics
} // namespace hydra

#endif
