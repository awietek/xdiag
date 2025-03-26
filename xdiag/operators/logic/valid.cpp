#include "valid.hpp"

#include <set>
#include <string>

#include <xdiag/operators/logic/types.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

void check_valid(Op const &op) try {
  std::string type = op.type();
  if (is_known_type(type)) {
    if ((type == "Sz") || (type == "S+") || (type == "S-") ||
        (type == "Cdagup") || (type == "Cup") || (type == "Cdagdn") ||
        (type == "Cdn") || (type == "Ntot") || (type == "Nup") ||
        (type == "Ndn") || (type == "Nupdn")) {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_nsites(op, 1);
    } else if ((type == "Hop") || (type == "Hopup") || (type == "Hopdn") ||
               (type == "SzSz") || (type == "SdotS") || (type == "Exchange") ||
               (type == "NtotNtot") || (type == "NupdnNupdn") ||
               (type == "tJSzSz") || (type == "tJSdotS")) {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_nsites(op, 2);
    } else if (type == "ScalarChirality") {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_nsites(op, 3);
      must_have_disjoint_sites(op);
    } else if ((type == "HubbardU") || (type == "Id")) {
      must_not_have_matrix(op);
      must_not_have_sites(op);
    } else if (type == "Matrix") {
      must_have_matrix(op);
      must_have_disjoint_sites(op);
    } else {
      XDIAG_THROW(
          "Logic error checking validity of Op (this is a bug, please report)");
    }
  } else {
    XDIAG_THROW(fmt::format("Unknown Op type: \"{}\", got Op:\n{}\nXDiag "
                            "recognizes the following Op types:\n{}",
                            type, to_string(op), known_types_string()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_valid(OpSum const &ops) try {
  for (auto [cpl, op] : ops) {
    check_valid(op);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_valid(Op const &op, int64_t nsites) try {
  check_valid(op);
  if (op.hassites()) {
    must_have_sites_in_range(op, 0, nsites);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_valid(OpSum const &ops, int64_t nsites) try {
  for (auto [cpl, op] : ops) {
    check_valid(op, nsites);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_sites(Op const &op) try {
  if (!op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_not_have_sites(Op const &op) try {
  if (op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" cannot have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_nsites(Op const &op, int64_t n) try {
  if (op.hassites()) {
    if (op.sites().size() != n) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have exactly {} sites defined, got Op:\n{}",
          op.type(), n, to_string(op)));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_disjoint_sites(Op const &op) try {
  if (op.hassites()) {
    auto const &sites = op.sites();
    auto set = std::set<int64_t>(sites.begin(), sites.end());
    if (set.size() != sites.size()) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have strictly disjoint sites, got Op:\n{}",
          op.type(), to_string(op)));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
void must_have_sites_in_range(Op const &op, int64_t l, int64_t u) try {
  if (op.hassites()) {
    for (auto s : op.sites()) {
      if ((s < 0) || (s >= u)) {
        XDIAG_THROW(fmt::format(
            "Op of type \"{}\" has site with index {}, but the "
            "indices must lie in the interval [{}, {}). Got Op:\n{}",
            op.type(), s, l, u, to_string(op)));
      }
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_matrix(Op const &op) try {
  if (!op.hasmatrix()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have a matrix defined, got Op:\n{}",
                    op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_not_have_matrix(Op const &op) try {
  if (op.hasmatrix()) {
    XDIAG_THROW(fmt::format(
        "Op of type \"{}\" cannot have a matrix defined, got Op:\n{}",
        op.type(), to_string(op)));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
