// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {

State apply(OpSum const &ops, State const &v) try {
  // invalid v
  if (!isvalid(v)) {
    return State();
  }

  // ops is zero
  if (isapprox(ops, OpSum())) {
    return State();
  }

  auto blockr = block(ops, v.block());
  bool real = isreal(ops) && isreal(v);
  auto w = State(blockr, real, v.ncols());
  apply(ops, v, w);
  return w;
}
XDIAG_CATCH

State apply(Op const &op, State const &v) try { return apply(OpSum(op), v); }
XDIAG_CATCH

void apply(OpSum const &ops, State const &v, State &w) try {

  if (!isvalid(v)) {
    w = State();
  } else if (isapprox(ops, OpSum())) { // ops is zero
    w = State();
  } else {

    if (!blocks_match(ops, v.block(), w.block())) {
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

void apply(Op const &op, State const &v, State &w) try {
  apply(OpSum(op), v, w);
}
XDIAG_CATCH

template <typename mat_t>
void apply(OpSum const &ops, Block const &block_in, mat_t const &mat_in,
           Block const &block_out, mat_t &mat_out) try {
  std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
          [&](tJ const &b1, tJ const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
          [&](Electron const &b1, Electron const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &b1, SpinhalfDistributed const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
          [&](tJDistributed const &b1, tJDistributed const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
          [&](ElectronDistributed const &b1, ElectronDistributed const &b2) {
            apply(ops, b1, mat_in, b2, mat_out);
          },
#endif
          [](auto const &, auto const &) {
            XDIAG_THROW(fmt::format("Invalid combination of Block types"));
          }},
      block_in, block_out);
}
XDIAG_CATCH

template void apply(OpSum const &, Block const &, arma::vec const &,
                    Block const &, arma::vec &);
template void apply(OpSum const &, Block const &, arma::cx_vec const &,
                    Block const &, arma::cx_vec &);
template void apply(OpSum const &, Block const &, arma::mat const &,
                    Block const &, arma::mat &);
template void apply(OpSum const &, Block const &, arma::cx_mat const &,
                    Block const &, arma::cx_mat &);

template <typename coeff_t>
static inline void fill_apply(coeff_t const *vec_in, coeff_t *vec_out,
                              int64_t idx_in, int64_t idx_out, coeff_t val) {
  // Atomic update to avoid multiple threads writing to the same address
#ifdef _OPENMP
  if constexpr (isreal<coeff_t>()) {
    coeff_t x = val * vec_in[idx_in];
#pragma omp atomic update
    vec_out[idx_out] += x;
  } else {
    complex x = val * vec_in[idx_in];
    double *r = &reinterpret_cast<double (&)[2]>(vec_out[idx_out])[0];
    double *i = &reinterpret_cast<double (&)[2]>(vec_out[idx_out])[1];
#pragma omp atomic update
    *r += x.real();
#pragma omp atomic update
    *i += x.imag();
  }
#else
  vec_out[idx_out] += val * vec_in[idx_in];
#endif
}

template <typename coeff_t>
static inline void fill_apply(arma::Col<coeff_t> const &vec_in,
                              arma::Col<coeff_t> &vec_out, int64_t idx_in,
                              int64_t idx_out, coeff_t val) {
  fill_apply(vec_in.memptr(), vec_out.memptr(), idx_in, idx_out, val);
}

template <typename coeff_t>
static inline void fill_apply(arma::Mat<coeff_t> const &mat_in,
                              arma::Mat<coeff_t> &mat_out, int64_t idx_in,
                              int64_t idx_out, coeff_t val) {
  // for each column call the usual fill_apply.
  for (int i = 0; i < mat_in.n_cols; i++) {
    fill_apply(mat_in.colptr(i), mat_out.colptr(i), idx_in, idx_out, val);
  }
}

template <typename mat_t, typename block_t>
void apply(OpSum const &ops, block_t const &block_in, mat_t const &mat_in,
           block_t const &block_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;

  check_valid(ops, block_in.nsites());
  mat_out.zeros();
  OpSum opsc = operators::compile<block_t>(ops);

  if constexpr (isdistributed<block_t>()) {
    algebra::apply_dispatch<coeff_t>(opsc, block_in, mat_in, block_out,
                                     mat_out);
  } else {
    // create fill method to apply the terms
    auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
      return fill_apply(mat_in, mat_out, idx_in, idx_out, val);
    };
#ifdef _OPENMP
    omp_set_schedule(omp_sched_guided, 0);
#endif
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  }
}
XDIAG_CATCH

template void apply(OpSum const &, Spinhalf const &, arma::vec const &,
                    Spinhalf const &, arma::vec &);
template void apply(OpSum const &, Spinhalf const &, arma::cx_vec const &,
                    Spinhalf const &, arma::cx_vec &);
template void apply(OpSum const &, Spinhalf const &, arma::mat const &,
                    Spinhalf const &, arma::mat &);
template void apply(OpSum const &, Spinhalf const &, arma::cx_mat const &,
                    Spinhalf const &, arma::cx_mat &);

template void apply(OpSum const &, tJ const &, arma::vec const &, tJ const &,
                    arma::vec &);
template void apply(OpSum const &, tJ const &, arma::cx_vec const &, tJ const &,
                    arma::cx_vec &);
template void apply(OpSum const &, tJ const &, arma::mat const &, tJ const &,
                    arma::mat &);
template void apply(OpSum const &, tJ const &, arma::cx_mat const &, tJ const &,
                    arma::cx_mat &);

template void apply(OpSum const &, Electron const &, arma::vec const &,
                    Electron const &, arma::vec &);
template void apply(OpSum const &, Electron const &, arma::cx_vec const &,
                    Electron const &, arma::cx_vec &);
template void apply(OpSum const &, Electron const &, arma::mat const &,
                    Electron const &, arma::mat &);
template void apply(OpSum const &, Electron const &, arma::cx_mat const &,
                    Electron const &, arma::cx_mat &);

#ifdef XDIAG_USE_MPI
template void apply(OpSum const &, SpinhalfDistributed const &,
                    arma::vec const &, SpinhalfDistributed const &,
                    arma::vec &);
template void apply(OpSum const &, SpinhalfDistributed const &,
                    arma::cx_vec const &, SpinhalfDistributed const &,
                    arma::cx_vec &);
template void apply(OpSum const &, SpinhalfDistributed const &,
                    arma::mat const &, SpinhalfDistributed const &,
                    arma::mat &);
template void apply(OpSum const &, SpinhalfDistributed const &,
                    arma::cx_mat const &, SpinhalfDistributed const &,
                    arma::cx_mat &);

template void apply(OpSum const &, tJDistributed const &, arma::vec const &,
                    tJDistributed const &, arma::vec &);
template void apply(OpSum const &, tJDistributed const &, arma::cx_vec const &,
                    tJDistributed const &, arma::cx_vec &);
template void apply(OpSum const &, tJDistributed const &, arma::mat const &,
                    tJDistributed const &, arma::mat &);
template void apply(OpSum const &, tJDistributed const &, arma::cx_mat const &,
                    tJDistributed const &, arma::cx_mat &);

template void apply(OpSum const &, ElectronDistributed const &,
                    arma::vec const &, ElectronDistributed const &,
                    arma::vec &);
template void apply(OpSum const &, ElectronDistributed const &,
                    arma::cx_vec const &, ElectronDistributed const &,
                    arma::cx_vec &);
template void apply(OpSum const &, ElectronDistributed const &,
                    arma::mat const &, ElectronDistributed const &,
                    arma::mat &);
template void apply(OpSum const &, ElectronDistributed const &,
                    arma::cx_mat const &, ElectronDistributed const &,
                    arma::cx_mat &);
#endif

} // namespace xdiag
