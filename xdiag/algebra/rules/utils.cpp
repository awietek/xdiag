// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "utils.hpp"

#include <vector>

namespace xdiag::operators {

OpSum replace_pair(Monomial const &mono, int64_t k, OpSum const &repl_pair) {
  std::vector<Op> pre(mono.ops().begin(), mono.ops().begin() + k);
  std::vector<Op> suf(mono.ops().begin() + k + 2, mono.ops().end());
  Monomial prefix(pre), suffix(suf);

  OpSum result;
  for (auto const &[c, m] : repl_pair) {
    result += OpSum(c, prefix * m * suffix);
  }
  return result;
}

OpSum swap_pair(Monomial const &mono, int64_t k, double sign) {
  std::vector<Op> pre(mono.ops().begin(), mono.ops().begin() + k);
  std::vector<Op> suf(mono.ops().begin() + k + 2, mono.ops().end());
  Monomial prefix(pre), suffix(suf);
  Monomial swapped({mono[k + 1], mono[k]});
  return sign * (prefix * swapped * suffix);
}

} // namespace xdiag::operators
