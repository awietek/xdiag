// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "expect.hpp"

#include <optional>
#include <variant>

#include <xdiag/algebra/symmetrize.hpp>
#include <xdiag/blocks/blocks.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/types.hpp>
#include <xdiag/states/apply.hpp>
#include <xdiag/states/dot.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation_set.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

namespace {

// Shared implementation of expect / expectC. coeff_t (double or complex)
// selects the real dot / complex dotC used to contract <state| sym(Op_i) |w>.
template <typename coeff_t>
arma::Col<coeff_t> expect_static(State const &state, std::string const &type) {
  int64_t n = state.nsites();
  if (n == 0) {
    return arma::Col<coeff_t>();
  }
  Block blk = state.block();

  // On a symmetry-adapted block a bare single-site Op has no well-defined
  // permutation irrep; symmetrizing over the site-permutation group (trivial,
  // all-ones characters) makes it group-invariant, so it acts within the block
  // and <state| sym(Op_i) |state> equals the per-site expectation.
  std::optional<PermutationGroup> group =
      std::visit([](auto const &b) { return b.irreps(); }, blk)
          .group("SitePermutation");
  auto site_op = [&](int64_t i) -> OpSum {
    return group ? symmetrize(Op(type, i), *group) : OpSum(Op(type, i));
  };

  // The operator's quantum-number shift is site-independent. If it changes the
  // block's quantum numbers, Op_i maps state into an orthogonal block and the
  // expectation is exactly zero on every site -- no computation needed.
  if (!blocks_match(site_op(0), blk, blk)) {
    return arma::Col<coeff_t>(n, arma::fill::zeros);
  }
  // Same block: reuse it (shallow copy) for the output state w, which the
  // in-place apply then fills for every site.
  State w(blk, state.isreal());

  arma::Col<coeff_t> result(n);
  for (int64_t i = 0; i < n; ++i) {
    apply(site_op(i), state, w);
    if constexpr (isreal<coeff_t>()) {
      result(i) = dot(state, w);
    } else {
      result(i) = dotC(state, w);
    }
  }
  return result;
}

} // namespace

arma::vec expect(State const &state, std::string type) try {
  if (!isreal(state) || !is_real_type(type)) {
    XDIAG_THROW("Cannot compute a real expectation value: the state or the "
                "operator \"" +
                type +
                "\" is complex, so the result may be complex. Consider using "
                "expectC instead.");
  }
  return expect_static<double>(state, type);
}
XDIAG_CATCH

arma::cx_vec expectC(State const &state, std::string type) try {
  return expect_static<complex>(state, type);
}
XDIAG_CATCH

} // namespace xdiag
