#include "order.hpp"

#include <algorithm>
#include <xdiag/operators/logic/types.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {

static int64_t permuted_index(int64_t idx, int64_t d,
                              std::vector<int64_t> const &perm) {
  std::vector<int64_t> perm_inv(perm.size());
  for (int64_t i = 0; i < (int64_t)perm.size(); ++i) {
    perm_inv[perm[i]] = i;
  }

  int64_t exp = 0;
  int64_t idx_permuted = 0;
  while (idx) {
    int64_t local = idx % d;
    idx_permuted += local * pow(d, perm_inv[exp++]);
    idx /= d;
  }

  return idx_permuted;
}

template <typename T>
static arma::Mat<T> permute_matrix(arma::Mat<T> const &mat,
                                   std::vector<int64_t> perm) try {
  int64_t m = mat.n_rows;
  int64_t n = mat.n_cols;
  if (m != n) {
    XDIAG_THROW(fmt::format("Matrix is not square, size: ({}, {})", m, n));
  }

  // Determine local dimension of matrix
  int64_t d = 0;
  int64_t nsites = perm.size();
  while (pow(d, nsites) < m) {
    ++d;
  }

  if (pow(d, nsites) != m) {
    XDIAG_THROW(fmt::format("Matrix dimensions are not of the form d^N, where "
                            "N denotes the number of sites (here {})",
                            nsites));
  }

  // Compute the permuted matrix
  arma::Mat<T> mat_permuted(m, n);
  for (int64_t i = 0; i < m; ++i) {
    int64_t ip = permuted_index(i, d, perm);
    for (int64_t j = 0; j < m; ++j) {
      int64_t jp = permuted_index(j, d, perm);
      mat_permuted(ip, jp) = mat(i, j);
    }
  }
  return mat_permuted;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::pair<Scalar, Op> order(Op const &op) try {
  return order(1.0, op);
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

  } else if (nsites_of_type(type) == 2) {
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
    if ((s1 < s2) && (s2 < s3)) {
      return {alpha, Op(type, {s1, s2, s3})};
    } else if ((s2 < s3) && (s3 < s1)) {
      return {alpha, Op(type, {s2, s3, s1})};
    } else if ((s3 < s1) && (s1 < s2)) {
      return {alpha, Op(type, {s3, s1, s2})};
    }
    // otherwise, we get a negative sign
    else if ((s1 < s3) && (s3 < s2)) {
      return {-alpha, Op(type, {s1, s3, s2})};
    } else if ((s3 < s2) && (s2 < s1)) {
      return {-alpha, Op(type, {s3, s2, s1})};
    } else { //  ((s2 < s1) || (s1 < s3))
      return {-alpha, Op(type, {s2, s1, s3})};
    }
  } else {
    return {alpha, op};
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

static bool less_sites(std::vector<int64_t> const &s1,
                       std::vector<int64_t> const &s2) try {
  int64_t n1 = s1.size();
  int64_t n2 = s2.size();
  if (n1 != n2) {
    return n1 < n2;
  } else { // Both number of sites are equal -> lexicographic site order
    return std::lexicographical_compare(s1.begin(), s1.end(), s2.begin(),
                                        s2.end());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T>
static bool less_matrix(arma::Mat<T> const &mat1, arma::Mat<T> const &mat2) {
  int64_t m1 = mat1.n_rows;
  int64_t n1 = mat1.n_cols;
  int64_t m2 = mat2.n_rows;
  int64_t n2 = mat2.n_cols;
  if ((m1 == m2) && (n1 == n2)) {

    // Elementwise lexicographic ordering
    if constexpr (isreal<T>()) {
      int64_t s = m1 * n1;
      const double *p1 = mat1.memptr();
      const double *p2 = mat2.memptr();
      return std::lexicographical_compare(p1, p1 + s, p2, p2 + s);
    } else {
      int64_t s = m1 * n1 * 2;
      const double *p1 = reinterpret_cast<const double *>(mat1.memptr());
      const double *p2 = reinterpret_cast<const double *>(mat2.memptr());
      return std::lexicographical_compare(p1, p1 + s, p2, p2 + s);
    }

  } else {
    return (m1 < m2) || ((m1 == m2) && (n1 < n2));
  }
}

bool less(Matrix const &mat1, Matrix const &mat2) try {
  if (mat1.isreal() && mat2.isreal()) {
    return less_matrix(mat1.as<arma::mat>(), mat2.as<arma::mat>());
  } else if (!mat1.isreal() && !mat2.isreal()) {
    return less_matrix(mat1.as<arma::cx_mat>(), mat2.as<arma::cx_mat>());
  } else if (mat1.isreal() && !mat2.isreal()) {
    return true;
  } else {
    return false;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool less(Op const &o1, Op const &o2) try {
  std::string type1 = o1.type();
  std::string type2 = o2.type();
  if (type1 != type2) {
    return type1 < type2;
  } else { // type equal

    if (o1.hassites() && o2.hassites()) {
      auto const &s1 = o1.sites();
      auto const &s2 = o2.sites();
      if (s1 != s2) {
        return less_sites(s1, s2);
      } else { // sites are equal
        if (o1.hasmatrix() && o2.hasmatrix()) {
          Matrix const &mat1 = o1.matrix();
          Matrix const &mat2 = o2.matrix();
          return less(mat1, mat2);
        } else if (!o1.hasmatrix() && o2.hasmatrix()) {
          return true;
        } else {
          return false;
        }
      }
    } else if (!o1.hassites() && o2.hassites()) {
      return true;
    } else {
      return false;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool less(Scalar const &s1, Scalar const &s2) {
  if (s1.isreal() && s2.isreal()) {
    return s1.real() < s2.real();
  } else if (!s1.isreal() && !s2.isreal()) {
    return (s1.real() < s2.real()) ||
           ((s1.real() == s2.real()) && (s1.imag() < s2.imag()));
  } else if (s1.isreal() && !s2.isreal()) {
    return true;
  } else {
    return false;
  }
}

OpSum order(OpSum const &ops) try {
  using pair_t = std::pair<Scalar, Op>;
  std::vector<pair_t> terms;
  for (auto [a, op] : ops.plain()) {
    terms.push_back(order(a.scalar(), op));
  }
  std::sort(terms.begin(), terms.end(), [](pair_t const &p1, pair_t const &p2) {
    Scalar const &c1 = p1.first;
    Scalar const &c2 = p2.first;
    Op const &o1 = p1.second;
    Op const &o2 = p2.second;
    return less(o1, o2) || ((o1 == o2) && less(c1, c2));
  });
  OpSum ops_ordered;
  for (auto const &[s, o] : terms) {
    ops_ordered += s * o;
  }
  return ops_ordered;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
