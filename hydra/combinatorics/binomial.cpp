#include "binomial.h"
#include <array>

namespace hydra::combinatorics {

// Create a table for 32 x 32 binomial coefficients
static constexpr int max_d = 5;
static constexpr int max_nk = (1 << max_d);
static constexpr int max_nk2 = max_nk * max_nk;

static constexpr std::array<idx_t, max_nk2> get_binomial_table() {
  std::array<idx_t, max_nk2> table{0};
  for (int n = 0; n < max_nk; ++n)
    for (int k = 0; k < max_nk; ++k)
      table[(n << max_d) | k] = binom(n, k);
  return table;
}

static constexpr auto binomial_table = get_binomial_table();

idx_t binomial(int n, int k) {

  // compute binomial if n,k >= 32
  if ((n >> max_d) || (k >> max_d))
    return binom(n, k);

  // otherwise look it up
  else
    return binomial_table[(n << max_d) | k];
}

} // namespace hydra::combinatorics
