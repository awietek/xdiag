#include "apply.hpp"

#include <variant>

#include <xdiag/algebra/fill.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

#include <xdiag/basis/electron/apply/dispatch_apply.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch_apply.hpp>
#include <xdiag/basis/tj/apply/dispatch_apply.hpp>
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/apply/dispatch_apply.hpp>
#include <xdiag/basis/tj_distributed/apply/dispatch_apply.hpp>
#endif

namespace xdiag {

State apply(OpSum const &ops, State const &v) try {
  auto blockr = block(ops, v.block());
  bool real = isreal(ops) && isreal(v);
  auto w = State(blockr, real, v.ncols());
  apply(ops, v, w);
  return w;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State apply(Op const &op, State const &v) try {
  return apply(OpSum(op), v);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

void apply(OpSum const &ops, State const &v, State &w) try {
  if (!blocks_match(ops, v.block(), w.block())) {
    XDIAG_THROW("Cannot apply OpSum to State. The resulting state is not in "
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
    XDIAG_THROW("Applying a OpSum to a state with multiple "
                "columns generically not yet implemented "
                "(are the States of the same size?)");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

void apply(Op const &op, State const &v, State &w) try {
  apply(OpSum(op), v, w);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
#endif
          [](auto const &, auto const &) {
            XDIAG_THROW(fmt::format("Invalid combination of Block types"));
          }},
      block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply(OpSum const &, Block const &, arma::vec const &,
                    Block const &, arma::vec &);
template void apply(OpSum const &, Block const &, arma::cx_vec const &,
                    Block const &, arma::cx_vec &);
template void apply(OpSum const &, Block const &, arma::mat const &,
                    Block const &, arma::mat &);
template void apply(OpSum const &, Block const &, arma::cx_mat const &,
                    Block const &, arma::cx_mat &);

template <typename mat_t, typename block_t>
void apply(OpSum const &ops, block_t const &block_in, mat_t const &mat_in,
           block_t const &block_out, mat_t &mat_out) try {
  check_valid(ops, block_in.nsites());
  mat_out.zeros();
  OpSum opsc = operators::compile<block_t>(ops);
  basis::dispatch_apply(opsc, block_in, mat_in, block_out, mat_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

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
#endif

} // namespace xdiag
