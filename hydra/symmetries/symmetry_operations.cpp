#include "symmetry_operations.h"

#include <algorithm>
#include <vector>

#include <hydra/utils/bitops.h>

namespace hydra {
namespace detail {

bool is_valid_permutation(std::vector<int> const &permutation) {
  for (int i = 0; i < (int)permutation.size(); ++i) {
    if (std::find(permutation.begin(), permutation.end(), i) ==
        permutation.end())
      return false;
  }
  return true;
}

template <class bit_t>
bit_t apply_permutation(bit_t state, int n_sites, const int *permutation) {
  bit_t tstate = 0;
  for (int site = 0; site < n_sites; ++site) {
    tstate |= ((state >> site) & 1) << permutation[site];
  }
  return tstate;
}

template uint16 apply_permutation<uint16>(uint16, int, const int *);
template uint32 apply_permutation<uint32>(uint32, int, const int *);
template uint64 apply_permutation<uint64>(uint64, int, const int *);

template <class bit_t>
double fermi_sign(bit_t state, int n_sites, const int *permutation) {
  std::vector<int> sort(n_sites); // Bad for performance
  int sum = 0;
  int nf, k;

  for (nf = 0, k = 0; k < n_sites; ++k)
    if (utils::gbit(state, k)) {
      sort[nf] = permutation[k];
      ++nf;
    }

  int old_sum;
  while (true) {
    old_sum = sum;
    for (k = 0; k < (nf - 1); ++k)
      if (sort[k + 1] < sort[k]) {
        sum++;
        std::swap(sort[k + 1], sort[k]);
      }
    if (old_sum == sum)
      break;
  }
  return (sum % 2 ? -1 : 1);
}

template double fermi_sign<uint16>(uint16, int, const int *);
template double fermi_sign<uint32>(uint32, int, const int *);
template double fermi_sign<uint64>(uint64, int, const int *);

} // namespace detail
} // namespace hydra
