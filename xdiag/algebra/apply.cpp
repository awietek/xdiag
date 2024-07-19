#include "apply.hpp"

#include <variant>
#include <xdiag/blocks/electron/apply.hpp>
#include <xdiag/blocks/spinhalf/apply.hpp>
#include <xdiag/blocks/tj/apply.hpp>
#include <xdiag/common.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/blocks/spinhalf_distributed/apply.hpp>
#include <xdiag/blocks/tj_distributed/apply.hpp>
#endif

namespace xdiag {

template <typename coeff_t>
void apply(OpSum const &ops, Block const &block_in,
           arma::Col<coeff_t> const &vec_in, Block const &block_out,
           arma::Col<coeff_t> &vec_out, double zero_precision) try {
  std::visit(
      [&](auto &&block_in, auto &&block_out) {
        apply(ops, block_in, vec_in, block_out, vec_out, zero_precision);
      },
      block_in, block_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply(OpSum const &, Block const &, arma::vec const &,
                    Block const &, arma::vec &, double);
template void apply(OpSum const &, Block const &, arma::cx_vec const &,
                    Block const &, arma::cx_vec &, double);

void apply(OpSum const &ops, State const &v, State &w,
           double zero_precision) try {
  if ((v.n_cols() == 1) && (w.n_cols() == 1)) {
    if (ops.isreal()) {
      if (v.isreal() && w.isreal()) {
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w.vector(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
      } else if (v.isreal() && !w.isreal()) {
        auto w2 = State(w.block(), true);
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w2.vector(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
        w = w2;
      } else if (!v.isreal() && w.isreal()) {
        w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
      } else if (!v.isreal() && !w.isreal()) {
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
      }
    } else {
      if (v.isreal()) {
        auto v2 = v;
        v2.make_complex();
        w.make_complex();
        arma::cx_vec vvec = v2.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
      } else {
        w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, zero_precision);
      }
    }
  } else {
    XDIAG_THROW("Applying a OpSum to a state with multiple "
                "columns not yet implemented");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
