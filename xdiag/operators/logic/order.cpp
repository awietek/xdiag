#include "order.hpp"

#include <algorithm>
#include <xdiag/operators/logic/types.hpp>

namespace xdiag {

template <typename T>
static arma::Mat<T> permute_matrix(arma::Mat<T> const &mat,
                                   std::vector<int64_t> perm) try {
  int64_t m = mat.n_rows;
  int64_t n = mat.n_cols;
  if (m != n) {
    XDIAG_THROW(fmt::format("Matrix is not square, size: ({}, {})", m, n));
  }

  int64_t n_sites = perm.size();

  // Determine local dimension of matrix
  int64_t d = 0;
  while (pow(d, n_sites) < m) {
    ++d;
  }

  if (pow(d, n_sites) != m) {
    XDIAG_THROW(fmt::format("Matrix dimensions are not of the form d^N, where "
                            "N denotes the number of sites (here {})",
                            n_sites));
  }

  for (int64_t i = 0; i < m; ++i) {
    for (int64_t j = 0; j < m; ++j) {
      std::vector<int64_t> inds(n_sites);


    }
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::pair<Scalar, Op> order(Scalar const &alpha, Op const &op) try {
  check_valid(op);

  std::string type = op.type();

  if (type == "Matrix") {
    std::vector<int64_t> sites = op.sites();
    std::vector<int64_t> sites_sorted(sites);
    std::sort(sites_sorted.begin(), sites_sorted.end());
    std::vector<int64_t> perm(sites.size());
    std::iota(perm.begin(), perm.end(), 0);
    std::sort(perm.begin(), perm.end(),
              [&](int64_t s1, int64_t s2) { return (sites[s1] < sites[s2]); });

    Matrix mat = op.matrix();
    if (mat.isreal()) {
      arma::mat matp = permute_matrix(mat.as<arma::mat>(), perm);
      return {alpha, Op(type, sites_sorted, matp)};
    } else {
      arma::cx_mat matp = permute_matrix(mat.as<arma::cx_mat>(), perm);
      return {alpha, Op(type, sites_sorted, matp)};
    }

  } else if (n_sites_of_type(type) == 2) {
    int64_t s1 = op[0];
    int64_t s2 = op[1];

    // No reordering required, just return as is
    if (s1 < s2) {
      return {alpha, op};
    }

    // Reordering required, some operators require conjugating the scalar
    else {
      if ((type == "Exchange") || (type == "Hop") || (type == "Hopup") ||
          (type == "Hopdn")) {
        return {alpha.conj(), Op(type, {s2, s1})};
      } else {
        return {alpha, Op(type, {s2, s1})};
      }
    }

  } else if (type == "ScalarChirality") {
    int64_t s1 = op[0];
    int64_t s2 = op[1];
    int64_t s3 = op[2];

    // Even permutations -> return as is
    if (((s1 < s2) && (s2 < s3)) || ((s2 < s3) && (s3 < s1)) ||
        ((s3 < s1) && (s1 < s2))) {
      return {alpha, op};
    }

    // otherwise, we get a negative sign
    else if ((s1 < s3) || (s3 < s2)) {
      return {-alpha, Op(type, {s1, s3, s2})};
    } else if ((s3 < s2) || (s2 < s1)) {
      return {-alpha, Op(type, {s3, s2, s1})};
    } else if ((s2 < s1) || (s1 < s3)) {
      return {-alpha, Op(type, {s2, s3, s1})};
    }
  } else {
    return {alpha, op};
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool less(Op const &o1, Op const &o2);
OpSum order(OpSum const &ops);

} // namespace xdiag
