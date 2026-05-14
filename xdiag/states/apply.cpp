// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/armadillo.hpp>
#include <xdiag/matrices/apply.hpp>
#include <xdiag/algebra/isapprox.hpp>
#include <xdiag/operators/qns/block.hpp>
#include <xdiag/operators/qns/blocks_match.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag {

template <typename op_t> State apply(op_t const &ops, State const &v) try {
  // invalid v
  if (!isvalid(v)) {
    return State();
  }

  // ops is zero
  if (isapprox(OpSum(ops), OpSum())) {
    return State();
  }

  auto blockr = block(OpSum(ops), v.block());
  bool real = isreal(ops) && isreal(v);
  auto w = State(blockr, real, v.ncols());
  apply(ops, v, w);
  return w;
}
XDIAG_CATCH

template <typename op_t>
void apply(op_t const &ops, State const &v, State &w) try {

  if (!isvalid(v)) {
    w = State();
  } else if (isapprox(OpSum(ops), OpSum())) { // ops is zero
    w = State();
  } else {

    if (!blocks_match(OpSum(ops), v.block(), w.block())) {
      XDIAG_THROW(
          "Cannot apply OpSum to State. The resulting state is not in "
          "the correct symmetry sector. Please check the quantum numbers "
          "of the output state w.");
    }

    if ((v.ncols() == 1) && (w.ncols() == 1)) {
      if (isreal(ops)) {
        if (isreal(v) && isreal(w)) {
          arma::vec vvec = v.vector(0, false);
          arma::vec wvec = w.vector(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
        } else if (isreal(v) && !isreal(w)) {
          auto w2 = State(w.block(), true);
          arma::vec vvec = v.vector(0, false);
          arma::vec wvec = w2.vector(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
          w = w2;
        } else if (!isreal(v) && isreal(w)) {
          w.make_complex();
          arma::cx_vec vvec = v.vectorC(0, false);
          arma::cx_vec wvec = w.vectorC(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
        } else if (!isreal(v) && !isreal(w)) {
          arma::cx_vec vvec = v.vectorC(0, false);
          arma::cx_vec wvec = w.vectorC(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
        }
      } else {
        if (isreal(v)) {
          auto v2 = v;
          v2.make_complex();
          w.make_complex();
          arma::cx_vec vvec = v2.vectorC(0, false);
          arma::cx_vec wvec = w.vectorC(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
        } else {
          w.make_complex();
          arma::cx_vec vvec = v.vectorC(0, false);
          arma::cx_vec wvec = w.vectorC(0, false);
          apply(ops, v.block(), vvec, w.block(), wvec);
        }
      }
    } else if (v.ncols() == w.ncols()) {
      if (isreal(ops)) {
        if (isreal(v) && isreal(w)) {
          arma::mat vmat = v.matrix(false);
          arma::mat wmat = w.matrix(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
        } else if (isreal(v) && !isreal(w)) {
          auto w2 = State(w.block(), true, w.ncols());
          arma::mat vmat = v.matrix(false);
          arma::mat wmat = w2.matrix(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
          w = w2;
        } else if (!isreal(v) && isreal(w)) {
          w.make_complex();
          arma::cx_mat vmat = v.matrixC(false);
          arma::cx_mat wmat = w.matrixC(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
        } else if (!isreal(v) && !isreal(w)) {
          arma::cx_mat vmat = v.matrixC(false);
          arma::cx_mat wmat = w.matrixC(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
        }
      } else {
        if (isreal(v)) {
          auto v2 = v;
          v2.make_complex();
          w.make_complex();
          arma::cx_mat vmat = v2.matrixC(false);
          arma::cx_mat wmat = w.matrixC(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
        } else {
          w.make_complex();
          arma::cx_mat vmat = v.matrixC(false);
          arma::cx_mat wmat = w.matrixC(false);
          apply(ops, v.block(), vmat, w.block(), wmat);
        }
      }
    } else {
      XDIAG_THROW(
          fmt::format("Applying an OpSum to a state with multiple columns "
                      "results in a state with an equal number of columns. "
                      "However, the input state number of columns ({}) differs "
                      "from the output state number of columns ({}).",
                      v.ncols(), w.ncols()));
    }
  }
}
XDIAG_CATCH

template State apply(Op const &, State const &);
template State apply(OpSum const &, State const &);
template void apply(Op const &, State const &, State &);
template void apply(OpSum const &, State const &, State &);

} // namespace xdiag
