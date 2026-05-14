#include "nup_ndn.hpp"

#include <algorithm>
#include <set>
#include <string>

#include <xdiag/armadillo.hpp>
#include <xdiag/bits/popcount.hpp>
#include <xdiag/algebra/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::operators {

template <typename T>
static std::optional<int64_t> nup_matrix(arma::Mat<T> const &mat,
                                         double precision = 1e-12) {
  std::set<int64_t> diffs;
  int64_t diff;
  for (uint64_t i = 0; i < mat.n_rows; ++i) {
    for (uint64_t j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        diff = bits::popcount(i) -
               bits::popcount(j); // needs generalization for d!=2
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
}
XDIAG_CATCH

std::optional<int64_t> nup(Monomial const &mono) try {
  int64_t total_nup = 0;
  for (auto op : mono) {
    auto n = nup(op);
    if (n) {
      total_nup += *n;
    } else {
      return std::nullopt;
    }
  }
  return total_nup;
}
XDIAG_CATCH

std::optional<int64_t> nup(OpSum const &ops) try {
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> nups;
    for (auto [cpl, mono] : ops) {
      nups.push_back(nup(mono));
    }

    // Check if all elements are the same
    if (std::adjacent_find(nups.begin(), nups.end(), std::not_equal_to<>()) ==
        nups.end()) {
      return nups[0];
    } else {
      return std::nullopt;
    }
  }
}
XDIAG_CATCH

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
}
XDIAG_CATCH

std::optional<int64_t> ndn(Monomial const &mono) try {
  int64_t total_ndn = 0;
  for (auto op : mono) {
    auto n = ndn(op);
    if (n) {
      total_ndn += *n;
    } else {
      return std::nullopt;
    }
  }
  return total_ndn;
}
XDIAG_CATCH

std::optional<int64_t> ndn(OpSum const &ops) try {
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> ndns;
    for (auto [cpl, mono] : ops) {
      ndns.push_back(ndn(mono));
    }

    // Check if all elements are the same
    if (std::adjacent_find(ndns.begin(), ndns.end(), std::not_equal_to<>()) ==
        ndns.end()) {
      return ndns[0];
    } else {
      return std::nullopt;
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::operators
