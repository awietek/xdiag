// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "correlation_matrix.hpp"

#include <optional>
#include <variant>
#include <vector>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/scalar.hpp>
#include <xdiag/operators/hc.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/types.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag {

namespace {

// Shared implementation of correlation_matrix / correlation_matrixC. coeff_t
// (double or complex) selects the real dot / complex dotC.
template <typename coeff_t>
arma::Mat<coeff_t> correlation_matrix_static(State const &state,
                                             std::string const &type1,
                                             std::string const &type2) {
  int64_t n = state.nsites();
  if (n == 0) {
    return arma::Mat<coeff_t>();
  }
  Block blk = state.block();

  // On a symmetry-adapted block a two-site product has no well-defined
  // permutation irrep; symmetrizing over the site-permutation group (trivial,
  // all-ones characters) makes it group-invariant so it acts within the block.
  std::optional<PermutationGroup> group =
      std::visit([](auto const &b) { return b.irreps(); }, blk)
          .group("SitePermutation");
  auto pair_op = [&](int64_t i, int64_t j) -> OpSum {
    OpSum bare = OpSum(Monomial{Op(type1, i), Op(type2, j)});
    return group ? symmetrize(bare, *group) : bare;
  };

  // The product's quantum-number shift is site-independent. If it changes the
  // block's quantum numbers, every C(i,j) maps state into an orthogonal block
  // and is exactly zero -- no computation needed.
  if (!blocks_match(pair_op(0, 0), blk, blk)) {
    return arma::Mat<coeff_t>(n, n, arma::fill::zeros);
  }

#ifdef XDIAG_DISTRIBUTED
  // Distributed blocks only apply single-operator terms, so the two-site
  // product cannot be applied as one operator. Evaluate
  //   C(i,j) = <state| Op(type1,i) Op(type2,j) |state>
  // by two sequential single-site applies instead (distributed blocks carry no
  // permutation symmetry, so no symmetrization is needed here).
  if (isdistributed(blk)) {
    arma::Mat<coeff_t> result(n, n);
    // The second apply writes into w, which is pre-allocated on `blk` (same
    // basis instance as state). The product preserves the block, so this is the
    // correct output block, and dot/dotC -- which require an identical block --
    // accept it. (apply zeroes its output, so w can be reused.)
    State w(blk, state.isreal());
    for (int64_t i = 0; i < n; ++i) {
      for (int64_t j = 0; j < n; ++j) {
        apply(Op(type1, i), apply(Op(type2, j), state), w);
        if constexpr (isreal<coeff_t>()) {
          result(i, j) = dot(state, w);
        } else {
          result(i, j) = dotC(state, w);
        }
      }
    }
    return result;
  }
#endif

  // Many of the n*n operators coincide up to a scalar prefactor -- by site
  // permutation symmetry, because the two operators commute, or as a hermitian
  // conjugate of one already seen. Group them: each (i,j) is assigned a
  // representative operator, a prefactor lambda, and a conjugate flag, so that
  //   op(i,j) = lambda * rep        -> value = lambda * <rep>,            or
  //   op(i,j) = lambda * rep^dagger -> value = lambda * conj(<rep>).
  // Only the distinct representatives are actually applied to the state.
  algebra::Algebra algebra = algebra::symmetry_algebra(blk);
  double const tol = 1e-12;

  std::vector<OpSum> reps;               // distinct representative operators
  std::vector<OpSum> reps_hc;            // their hermitian conjugates
  std::vector<int64_t> rep_of(n * n);    // representative index for each (i,j)
  std::vector<coeff_t> lambda_of(n * n); // op(i,j) = lambda * rep [^dagger]
  std::vector<char> conj_of(n * n, 0);   // 1 if op(i,j) matches rep^dagger

  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      OpSum op = pair_op(i, j);
      int64_t idx = i * n + j;
      int64_t found = -1;
      for (int64_t r = 0; r < (int64_t)reps.size(); ++r) {
        // op = lambda * rep  ->  <op> = lambda * <rep>
        if (std::optional<Scalar> mult =
                isapprox_multiple(op, reps[r], algebra, tol, tol)) {
          found = r;
          lambda_of[idx] = mult->as<coeff_t>();
          conj_of[idx] = 0;
          break;
        }
        // op = lambda * rep^dagger  ->  <op> = lambda * conj(<rep>)
        if (std::optional<Scalar> mult =
                isapprox_multiple(op, reps_hc[r], algebra, tol, tol)) {
          found = r;
          lambda_of[idx] = mult->as<coeff_t>();
          conj_of[idx] = 1;
          break;
        }
      }
      if (found < 0) {
        found = (int64_t)reps.size();
        reps.push_back(op);
        reps_hc.push_back(hc(op));
        lambda_of[idx] = coeff_t(1.0);
        conj_of[idx] = 0;
      }
      rep_of[idx] = found;
    }
  }

  // Compute <state| rep |state> once per distinct representative, reusing the
  // output state w (same block) across all applies.
  State w(blk, state.isreal());
  std::vector<coeff_t> rep_value(reps.size());
  for (int64_t r = 0; r < (int64_t)reps.size(); ++r) {
    apply(reps[r], state, w);
    if constexpr (isreal<coeff_t>()) {
      rep_value[r] = dot(state, w);
    } else {
      rep_value[r] = dotC(state, w);
    }
  }

  // Assemble the matrix from the representative values, conjugating those that
  // matched a representative's hermitian conjugate. (For a real result conj is
  // the identity, so the flag has no effect.)
  arma::Mat<coeff_t> result(n, n);
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      int64_t idx = i * n + j;
      coeff_t v = rep_value[rep_of[idx]];
      if constexpr (!isreal<coeff_t>()) {
        if (conj_of[idx]) {
          v = std::conj(v);
        }
      }
      result(i, j) = lambda_of[idx] * v;
    }
  }
  return result;
}

} // namespace

arma::mat correlation_matrix(State const &state, std::string type1,
                             std::string type2) try {
  if (!isreal(state) || !is_real_type(type1) || !is_real_type(type2)) {
    XDIAG_THROW("Cannot compute a real correlation matrix: the state or one of "
                "the operators \"" +
                type1 + "\", \"" + type2 +
                "\" is complex, so the result may be complex. Consider using "
                "correlation_matrixC instead.");
  }
  return correlation_matrix_static<double>(state, type1, type2);
}
XDIAG_CATCH

arma::cx_mat correlation_matrixC(State const &state, std::string type1,
                                 std::string type2) try {
  return correlation_matrix_static<complex>(state, type1, type2);
}
XDIAG_CATCH

} // namespace xdiag
