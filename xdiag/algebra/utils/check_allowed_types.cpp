// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "check_allowed_types.hpp"

#include <cstdint>
#include <string>

#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::algebra {

void check_allowed_types(OpSum const &ops, Algebra const &algebra) try {
  if (algebra.allowed_types.empty()) {
    return;
  }
  for (auto const &[coeff, mono] : ops) {
    for (auto const &op : mono) {
      if (algebra.allowed_types.count(op.type()) == 0) {
        std::string allowed;
        bool first = true;
        for (std::string const &t : algebra.allowed_types) {
          if (!first) {
            allowed += ", ";
          }
          allowed += t;
          first = false;
        }
        XDIAG_THROW(
            fmt::format("Op of type \"{}\" is not allowed in the {} algebra.\n"
                        "Allowed types: {}",
                        op.type(), algebra.name, allowed));
      }
    }
  }
}
XDIAG_CATCH

void check_sites_in_range(OpSum const &ops, Algebra const &algebra) try {
  int64_t nsites = algebra.nsites;
  for (auto const &[coeff, mono] : ops) {
    for (auto const &op : mono) {
      if (!op.hassites()) {
        continue;
      }
      for (int64_t s : op.sites()) {
        if ((s < 0) || (s >= nsites)) {
          XDIAG_THROW(fmt::format(
              "Op of type \"{}\" acts on site {}, which is out of range for "
              "the \"{}\" algebra defined on {} site(s) (valid range "
              "[0, {})).",
              op.type(), s, algebra.name, nsites, nsites));
        }
      }
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::algebra
