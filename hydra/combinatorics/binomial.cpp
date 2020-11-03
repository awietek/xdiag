#include "binomial.h"

namespace hydra {
namespace combinatorics {

int64 binomial(const int &n, const int &k) {
  if (k > n || k < 0)
    return 0;
  int64 res = 1;
  for (int i = 1; i <= k; i++)
    res = (res * (n - i + 1)) / i;
  return res;
}


} // namespace combinatorics
} // namespace hydra
