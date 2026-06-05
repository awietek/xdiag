// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "valid.hpp"

#include <set>

#include <xdiag/operators/types.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::operators {

void check_valid(Op const &op) try {
  std::string type = op.type();

  if (!is_known_type(type)) {
    XDIAG_THROW(fmt::format("Op type \"{}\" is unknown. Known types: {}", type,
                            known_types_string()));
  }

  auto const &info = info_of_type(type);

  if (info.site_required && !op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    type, to_string(op)));
  }

  if (op.hassites()) {
    if (info.nsites != undefined && (int64_t)op.sites().size() != info.nsites) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have exactly {} site(s), got {}. Op:\n{}",
          type, info.nsites, op.sites().size(), to_string(op)));
    }
    if (!info.allow_overlapping) {
      auto const &sites = op.sites();
      auto s = std::set<int64_t>(sites.begin(), sites.end());
      if (s.size() != sites.size()) {
        XDIAG_THROW(fmt::format(
            "Op of type \"{}\" must have distinct sites, got Op:\n{}", type,
            to_string(op)));
      }
    }
  }

  if (info.matrix_required && !op.hasmatrix()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have a matrix defined, got Op:\n{}",
                    type, to_string(op)));
  }

  if (!info.matrix_required && op.hasmatrix()) {
    XDIAG_THROW(fmt::format(
        "Op of type \"{}\" must not carry a matrix. Got Op:\n{}",
        type, to_string(op)));
  }
}
XDIAG_CATCH

void check_valid(Monomial const &mono) try {
  for (auto const &op : mono) {
    check_valid(op);
  }
}
XDIAG_CATCH

void check_valid(OpSum const &ops) try {
  for (auto const &[cpl, mono] : ops) {
    check_valid(mono);
  }
}
XDIAG_CATCH

void must_have_sites(Op const &op) try {
  if (!op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_not_have_sites(Op const &op) try {
  if (op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" cannot have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_have_nsites(Op const &op, int64_t n) try {
  if (op.hassites()) {
    if ((int64_t)op.sites().size() != n) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have exactly {} site(s), got Op:\n{}",
          op.type(), n, to_string(op)));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_have_disjoint_sites(Op const &op) try {
  if (op.hassites()) {
    auto const &sites = op.sites();
    auto s = std::set<int64_t>(sites.begin(), sites.end());
    if (s.size() != sites.size()) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have strictly disjoint sites, got Op:\n{}",
          op.type(), to_string(op)));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_have_sites_in_range(Op const &op, int64_t l, int64_t u) try {
  if (op.hassites()) {
    for (auto s : op.sites()) {
      if ((s < l) || (s >= u)) {
        XDIAG_THROW(fmt::format(
            "Op of type \"{}\" has site {}, but sites must lie in [{}, {}). "
            "Got Op:\n{}",
            op.type(), s, l, u, to_string(op)));
      }
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_have_matrix(Op const &op) try {
  if (!op.hasmatrix()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have a matrix defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
}
XDIAG_CATCH

void must_not_have_matrix(Op const &op) try {
  if (op.hasmatrix()) {
    XDIAG_THROW(fmt::format(
        "Op of type \"{}\" cannot have a matrix defined, got Op:\n{}",
        op.type(), to_string(op)));
  }
}
XDIAG_CATCH

} // namespace xdiag::operators
