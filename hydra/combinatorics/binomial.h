#pragma once

#include <hydra/common.h>

namespace hydra::combinatorics {

idx_t binomial(int n, int k);

constexpr idx_t binom(int n, int k) {
  if (k > n || k < 0)
    return 0;
  idx_t res = 1;
  for (int i = 1; i <= k; i++)
    res = (res * (n - i + 1)) / i;
  return res;
}

} // namespace hydra::combinatorics
