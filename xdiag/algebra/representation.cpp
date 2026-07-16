// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "representation.hpp"

#include <cmath>
#include <map>
#include <set>
#include <vector>

#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/permute.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/operators/valid.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::algebra {

// --- U(1) charge of an OpSum ------------------------------------------------

// Charge each elementary Op type carries under the named U(1) action. Types not
// listed carry charge 0.
static std::map<std::string, int64_t>
u1_charges(std::string const &action) try {
  if (action == "nup") {
    return {{"S+", 1}, {"Cdagup", 1}, {"S-", -1}, {"Cup", -1}};
  } else if (action == "ndn") {
    return {{"S-", 1}, {"Cdagdn", 1}, {"S+", -1}, {"Cdn", -1}};
  } else if (action == "number") {
    return {{"Cdag", 1}, {"C", -1}, {"Adag", 1}, {"A", -1}};
  } else {
    XDIAG_THROW("Unknown U(1) action \"" + action + "\"");
  }
}
XDIAG_CATCH

// Particle number of a local many-site index: each of the `nsites` sites holds
// an occupation in 0..d-1, encoded as a digit of `index` in base `d`. The
// charge is the sum of occupations. Reduces to popcount for d == 2.
static int64_t particle_number(uint64_t index, int64_t d, int64_t nsites) {
  int64_t n = 0;
  for (int64_t s = 0; s < nsites; ++s) {
    n += static_cast<int64_t>(index % static_cast<uint64_t>(d));
    index /= static_cast<uint64_t>(d);
  }
  return n;
}

// Charge of a "Matrix" Op for local dimension d (a property of the algebra):
// the matrix acts on `nsites` sites and must be d^nsites x d^nsites. A nonzero
// entry (i, j) shifts the particle number by
// particle_number(i) - particle_number(j); the Op has a well-defined charge
// only if all nonzero entries agree.
template <typename T>
static std::optional<int64_t> matrix_charge(arma::Mat<T> const &mat, int64_t d,
                                            int64_t nsites, double tol) try {
  int64_t dim = math::ipow(d, nsites);
  if (((int64_t)mat.n_rows != dim) || ((int64_t)mat.n_cols != dim)) {
    XDIAG_THROW("Matrix Op dimension does not match d^nsites for the algebra's "
                "local dimension d");
  }

  std::set<int64_t> diffs;
  int64_t diff = 0;
  for (uint64_t i = 0; i < mat.n_rows; ++i) {
    for (uint64_t j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > tol) {
        diff = particle_number(i, d, nsites) - particle_number(j, d, nsites);
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
XDIAG_CATCH

static std::optional<int64_t>
op_charge(Op const &op, std::map<std::string, int64_t> const &charges,
          int64_t d, double tol) try {
  operators::check_valid(op);
  std::string type = op.type();
  if (type == "Matrix") {
    Matrix mat = op.matrix();
    if (isreal(mat)) {
      return matrix_charge(mat.as<arma::mat>(), d, op.size(), tol);
    } else {
      return matrix_charge(mat.as<arma::cx_mat>(), d, op.size(), tol);
    }
  }
  auto it = charges.find(type);
  return (it != charges.end()) ? it->second : 0;
}
XDIAG_CATCH

static std::optional<int64_t> opsum_charge(OpSum const &ops,
                                           std::string const &action, int64_t d,
                                           double tol) try {
  std::map<std::string, int64_t> charges = u1_charges(action);
  std::optional<int64_t> total;
  for (auto const &[coeff, mono] : ops) {
    int64_t mono_charge = 0;
    for (Op const &op : mono) {
      std::optional<int64_t> q = op_charge(op, charges, d, tol);
      if (!q) {
        return std::nullopt;
      }
      mono_charge += *q;
    }
    if (!total) {
      total = mono_charge;
    } else if (*total != mono_charge) {
      return std::nullopt;
    }
  }
  return total ? total : std::optional<int64_t>(0);
}
XDIAG_CATCH

// --- representation() -------------------------------------------------------

// Determine the Representation that `ops` transforms under, with respect to the
// symmetry encoded by `irrep`. Only the symmetry of `irrep` is used: a
// permutation Representation supplies the PermutationGroup whose action is
// probed; a charge (U(1)) Representation supplies the action string via its
// type. The charge/characters carried by `irrep` itself are ignored and
// recomputed from `ops`. Returns std::nullopt if `ops` has no well-defined
// sector under that symmetry.
std::optional<Representation> representation(OpSum const &ops,
                                             Representation const &irrep,
                                             Algebra const &algebra,
                                             double tol) try {
  if (irrep.is_permutation()) {
    // `ops` transforms as a 1-D irrep iff every group element only rescales it.
    // The scalar lambda with permute(ops, g) == lambda * ops is the character
    // of g; if any element maps ops to a different operator there is no
    // well-defined sector.
    PermutationGroup group = irrep.group();
    int64_t n = group.size();
    arma::cx_vec chars(n);
    bool real = true;
    for (int64_t g = 0; g < n; ++g) {
      OpSum permuted = permute(ops, group[g]);
      std::optional<Scalar> lambda =
          isapprox_multiple(permuted, ops, algebra, tol, tol);
      if (!lambda) {
        return std::nullopt;
      }
      chars[g] = lambda->as<complex>();
      if (!lambda->isreal()) {
        real = false;
      }
    }
    if (real) {
      return Representation(group, arma::vec(arma::real(chars)));
    } else {
      return Representation(group, chars);
    }
  } else {
    // Charge (U(1)) symmetry: the type names the action ("nup", "ndn", "np").
    std::optional<int64_t> charge =
        opsum_charge(ops, irrep.type(), algebra.d, tol);
    if (!charge) {
      return std::nullopt;
    }
    return Representation(irrep.type(), *charge);
  }
}
XDIAG_CATCH

// Determine the Representation that `ops` transforms under for each symmetry in
// `irreps`. Symmetries under which `ops` has no well-defined sector are dropped
// from the result.
RepresentationSet representations(OpSum const &ops,
                                  RepresentationSet const &irreps,
                                  Algebra const &algebra, double tol) try {
  std::vector<Representation> result;
  for (Representation const &irrep : irreps) {
    std::optional<Representation> rep =
        representation(ops, irrep, algebra, tol);
    if (rep) {
      result.push_back(*rep);
    }
  }
  return RepresentationSet(result);
}
XDIAG_CATCH

} // namespace xdiag::algebra
