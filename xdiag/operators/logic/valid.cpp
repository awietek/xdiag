#include "valid.hpp"

#include <string>
#include <set>

#include <xdiag/operators/logic/types.hpp>

namespace xdiag {

void check_valid(Op const &op) try {
  std::string type = op.type();
  if (is_known_type(type)) {
    if ((type == "SZ") || (type == "S+") || (type == "S-") ||
        (type == "CDAGUP") || (type == "CUP") || (type == "CDAGDN") ||
        (type == "CDN") || (type == "NUMBER") || (type == "NUMBERUP") ||
        (type == "NUMBERDN")) {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_n_sites(op, 1);
    } else if ((type == "EXCHANGE") || (type == "ISING") || (type == "HOP") ||
               (type == "HOPUP") || (type == "HOPDN") || (type == "TJISING")) {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_n_sites(op, 2);
      must_have_disjoint_sites(op);
    } else if (type == "SCALARCHIRALITY") {
      must_not_have_matrix(op);
      must_have_sites(op);
      must_have_n_sites(op, 3);
      must_have_disjoint_sites(op);
    } else if (type == "HUBBARDU") {
      must_not_have_matrix(op);
      must_not_have_sites(op);
    } else if (type == "MATRIX") {
      must_have_matrix(op);
    } else {
      XDIAG_THROW("Logic error checking validity of Op (this is a bug)");
    }
  } else {
    XDIAG_THROW(fmt::format("Unknown Op type: \"{}\"", type));
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

void check_valid(Op const &op, int64_t n_sites) try {
  check_valid(op);
  must_have_sites_in_range(op, 0, n_sites);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void check_valid(OpSum const &ops, int64_t n_sites) try {
  for (auto [cpl, op] : ops) {
    check_valid(op, n_sites);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_sites(Op const &op) try {
  if (!op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined.", op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_not_have_sites(Op const &op) try {
  if (op.hassites()) {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" cannot have sites defined.", op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_n_sites(Op const &op, int64_t n) try {
  if (op.hassites()) {
    if (op.sites().size() != n) {
      XDIAG_THROW(fmt::format(
          "Op of type \"{}\" must have exactly {} sites defined.", n));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined.", op.type()));
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
          "Op of type \"{}\"must have strictly disjoint sites.", op.type()));
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined.", op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
void must_have_sites_in_range(Op const &op, int64_t l, int64_t u) try {
  if (op.hassites()) {
    for (auto s : op.sites()) {
      if ((s < 0) || (s >= u)) {
        XDIAG_THROW(
            fmt::format("Op of type \"{}\" has site with index {}, but the "
                        "indices must lie in the interval [{}, {}).",
                        op.type(), s, l, u));
      }
    }
  } else {
    XDIAG_THROW(
        fmt::format("Op of type \"{}\" must have sites defined.", op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_have_matrix(Op const &op) try {
  if (!op.hasmatrix()) {
    XDIAG_THROW(fmt::format("Op of type \"{}\" must have a matrix defined.",
                            op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

void must_not_have_matrix(Op const &op) try {
  if (op.hasmatrix()) {
    XDIAG_THROW(fmt::format("Op of type \"{}\" cannot have a matrix defined.",
                            op.type()));
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
