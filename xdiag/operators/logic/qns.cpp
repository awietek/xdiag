#include "qns.hpp"

#include <set>

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

template <typename T>
static std::optional<int64_t> nup_matrix(arma::Mat<T> const &mat,
                                         double precision = 1e-12) {
  std::set<int64_t> diffs;
  int64_t diff;
  for (uint64_t i = 0; i < mat.n_rows; ++i) {
    for (uint64_t j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        diff =
            bits::popcnt(i) - bits::popcnt(j); // needs generalization for d!=2
        diffs.insert(diff);
      }
    }
  }
  if (diffs.size() == 0) {
    return 0;
  } else if (diffs.size() == 1) {
    return diff;
  } else {
    return std::nullopt;
  }
}

std::optional<int64_t> nup(Op const &op) try {
  check_valid(op);
  std::string type = op.type();
  if ((type == "S+") || (type == "Cdagup")) {
    return 1;
  } else if ((type == "S-") || (type == "Cup")) {
    return -1;
  } else if (type == "Matrix") {
    Matrix mat = op.matrix();
    if (isreal(mat)) {
      return nup_matrix(mat.as<arma::mat>());
    } else {
      return nup_matrix(mat.as<arma::cx_mat>());
    }
  } else {
    return 0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> nup(OpSum const &ops) try {
  check_valid(ops);
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> nups;
    for (auto [cpl, op] : ops) {
      nups.push_back(nup(op));
    }

    // Check if all elements are the same
    if (std::adjacent_find(nups.begin(), nups.end(), std::not_equal_to<>()) ==
        nups.end()) {
      return nups[0];
    } else {
      return std::nullopt;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> ndn(Op const &op) try {
  check_valid(op);
  std::string type = op.type();
  if ((type == "S-") || (type == "Cdagdn")) {
    return 1;
  } else if ((type == "S+") || (type == "Cdn")) {
    return -1;
  } else {
    return 0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> ndn(OpSum const &ops) try {
  check_valid(ops);
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> ndns;
    for (auto [cpl, op] : ops) {
      ndns.push_back(ndn(op));
    }

    // Check if all elements are the same
    if (std::adjacent_find(ndns.begin(), ndns.end(), std::not_equal_to<>()) ==
        ndns.end()) {
      return ndns[0];
    } else {
      return std::nullopt;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<Representation>
representation(OpSum const &ops, PermutationGroup const &group) try {
  std::vector<complex> characters;
  std::vector<double> characters_real;
  bool real = true;
  for (auto const &perm : group) {
    OpSum opsp = permute(ops, perm);
    std::optional<Scalar> factor = isapprox_multiple(opsp, ops);
    if (factor) {
      characters.push_back((*factor).as<complex>());
      characters_real.push_back((*factor).real());
      if (!(*factor).isreal()) {
        real = false;
      }
    } else {
      return std::nullopt;
    }
  }
  if (real) {
    return Representation(group, characters_real);
  } else {
    return Representation(group, characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag
