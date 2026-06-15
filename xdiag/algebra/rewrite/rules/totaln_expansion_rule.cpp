// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "totaln_expansion_rule.hpp"

#include <optional>
#include <string>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

OpRule totaln_expansion_rule(int64_t nsites, std::string const &local_type) {
  return [nsites, local_type](Op const &op) -> std::optional<OpSum> {
    if (op.type() != "TotalN") {
      return std::nullopt;
    }
    // TotalN -> sum_i <local_type>{i}  (site-free: sums the local number
    // operator over all `nsites` sites; "N" for spinless fermion / boson,
    // "Ntot" for the spinful electron).
    OpSum r;
    for (int64_t i = 0; i < nsites; ++i) {
      r += Op(local_type, i);
    }
    return r;
  };
}

} // namespace xdiag::algebra
