#pragma once

#include <hydra/common.h>
#include <hydra/utils/bitops.h>

namespace hydra::combinatorics {

int64 binomial(int n, int k);

constexpr int64 binom(int n, int k) {
  if (k > n || k < 0)
    return 0;
  int64 res = 1;
  for (int i = 1; i <= k; i++)
    res = (res * (n - i + 1)) / i;
  return res;
}

} // namespace hydra::combinatorics
