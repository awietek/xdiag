// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "cyclic_group.hpp"

#include <vector>

#include <xdiag/math/numbers.hpp>
#include <xdiag/symmetries/permutation.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

PermutationGroup cyclic_group(int64_t n) try {

  if (n < 1) {
    XDIAG_THROW(fmt::format(
        "Invalid order of cyclic group: n={}, must be larger or equal to 1.",
        n));
  }

  std::vector<Permutation> permutation_array;
  for (int64_t sym = 0; sym < n; ++sym) {

    std::vector<int64_t> pv;
    for (int64_t site = 0; site < n; ++site) {
      int64_t newsite = (site + sym) % n;
      pv.push_back(newsite);
    }
    permutation_array.push_back(Permutation(pv));
  }
  return PermutationGroup(permutation_array);
}
XDIAG_CATCH

Representation cyclic_group_irrep(int64_t n, int64_t k) try {
  auto group = cyclic_group(n);
  if ((k < 0) || (k >= n)) {
    XDIAG_THROW(fmt::format(
        "Invalid irrep number: k={}, must be in interval [0, n) (given n={}).",
        k, n));
  }
  if (k == 0) { // real irrep, 0-momentum
    std::vector<double> chis(n, 1.0);
    return Representation(group, chis);
  } else if ((n % 2 == 0) && (k == n / 2)) { // real irrep, pi-momentum
    std::vector<double> chis(n);
    for (int64_t l = 0; l < n; ++l) {
      chis[l] = ((l % 2) == 0) ? 1.0 : -1.0;
    }
    return Representation(group, chis);
  } else { // complex irrep
    std::vector<complex> chis(n);
    for (int64_t l = 0; l < n; ++l) {
      chis[l] = complex{std::cos(2 * math::pi * l * k / n),
                        std::sin(2 * math::pi * l * k / n)};
    }
    return Representation(group, chis);
  }
}
XDIAG_CATCH

} // namespace xdiag
