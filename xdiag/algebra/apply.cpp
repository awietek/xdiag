#include "apply.hpp"

#include <variant>

#include <xdiag/algebra/fill.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/compiler.hpp>

#include <xdiag/basis/electron/apply/dispatch_apply.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch_apply.hpp>
#include <xdiag/basis/tj/apply/dispatch_apply.hpp>
#ifdef XDIAG_USE_MPI
#include <xdiag/basis/spinhalf_distributed/apply/dispatch_apply.hpp>
#include <xdiag/basis/tj_distributed/apply/dispatch_apply.hpp>
#endif

namespace xdiag {

void apply(Op const &op, State const &v, State &w, double precision) try {
  OpSum ops({op});
  apply(ops, v, w, precision);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

void apply(OpSum const &ops, State const &v, State &w, double precision) try {
  if ((v.n_cols() == 1) && (w.n_cols() == 1)) {
    if (ops.isreal()) {
      if (v.isreal() && w.isreal()) {
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w.vector(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
      } else if (v.isreal() && !w.isreal()) {
        auto w2 = State(w.block(), true);
        arma::vec vvec = v.vector(0, false);
        arma::vec wvec = w2.vector(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
        w = w2;
      } else if (!v.isreal() && w.isreal()) {
        w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
      } else if (!v.isreal() && !w.isreal()) {
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
      }
    } else {
      if (v.isreal()) {
        auto v2 = v;
        v2.make_complex();
        w.make_complex();
        arma::cx_vec vvec = v2.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
      } else {
        w.make_complex();
        arma::cx_vec vvec = v.vectorC(0, false);
        arma::cx_vec wvec = w.vectorC(0, false);
        apply(ops, v.block(), vvec, w.block(), wvec, precision);
      }
    }
  } else {
    XDIAG_THROW("Applying a OpSum to a state with multiple "
                "columns not yet implemented");
  }
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template <typename coeff_t>
void apply(OpSum const &ops, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_spinhalf(ops, n_sites, precision);
  vec_out.zeros();
  basis::spinhalf::dispatch_apply(opsc, block_in, vec_in, block_out, vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, Spinhalf const &,
                            arma::Col<double> const &, Spinhalf const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, Spinhalf const &,
                             arma::Col<complex> const &, Spinhalf const &,
                             arma::Col<complex> &, double);

template <typename coeff_t>
void apply(OpSum const &ops, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_tj(ops, n_sites, precision);
  vec_out.zeros();
  basis::tj::dispatch_apply(opsc, block_in, vec_in, block_out, vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, tJ const &,
                            arma::Col<double> const &, tJ const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, tJ const &,
                             arma::Col<complex> const &, tJ const &,
                             arma::Col<complex> &, double);

template <typename coeff_t>
void apply(OpSum const &ops, Electron const &block_in,
           arma::Col<coeff_t> const &vec_in, Electron const &block_out,
           arma::Col<coeff_t> &vec_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_electron(ops, n_sites, precision);
  vec_out.zeros();
  basis::electron::dispatch_apply(opsc, block_in, vec_in, block_out, vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, Electron const &,
                            arma::Col<double> const &, Electron const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, Electron const &,
                             arma::Col<complex> const &, Electron const &,
                             arma::Col<complex> &, double);

#ifdef XDIAG_USE_MPI

template <typename coeff_t>
void apply(OpSum const &ops, SpinhalfDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in,
           SpinhalfDistributed const &block_out, arma::Col<coeff_t> &vec_out,
           double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_spinhalf(ops, n_sites, precision);
  vec_out.zeros();
  basis::spinhalf_distributed::dispatch_apply(opsc, block_in, vec_in, block_out,
                                              vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, SpinhalfDistributed const &,
                            arma::Col<double> const &,
                            SpinhalfDistributed const &, arma::Col<double> &,
                            double);

template void apply<complex>(OpSum const &, SpinhalfDistributed const &,
                             arma::Col<complex> const &,
                             SpinhalfDistributed const &, arma::Col<complex> &,
                             double);

template <typename coeff_t>
void apply(OpSum const &ops, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out, double precision) try {
  int64_t n_sites = block_in.n_sites();
  OpSum opsc = operators::compile_tj(ops, n_sites, precision);
  vec_out.zeros();
  basis::tj_distributed::dispatch_apply(opsc, block_in, vec_in, block_out,
                                        vec_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template void apply<double>(OpSum const &, tJDistributed const &,
                            arma::Col<double> const &, tJDistributed const &,
                            arma::Col<double> &, double);

template void apply<complex>(OpSum const &, tJDistributed const &,
                             arma::Col<complex> const &, tJDistributed const &,
                             arma::Col<complex> &, double);

#endif

template <typename coeff_t>
void apply(OpSum const &ops, Block const &block_in,
           arma::Col<coeff_t> const &vec_in, Block const &block_out,
           arma::Col<coeff_t> &vec_out, double precision) try {
  std::visit(
      [&](auto &&block_in, auto &&block_out) {
        apply(ops, block_in, vec_in, block_out, vec_out, precision);
      },
      block_in, block_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply(OpSum const &, Block const &, arma::vec const &,
                    Block const &, arma::vec &, double);
template void apply(OpSum const &, Block const &, arma::cx_vec const &,
                    Block const &, arma::cx_vec &, double);

} // namespace xdiag
