// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// Dense matrix builder: fills a caller-allocated, Julia-owned column-major
// buffer of size size(out) x dim(in). The output block is inferred via
// block(ops, in) (so Julia never round-trips the Block variant). Templated on
// coeff_t, so instantiated here per block. Sparse builders live in sparse.cpp.

#include <julia/src/xdiagjl.hpp>
#include <xdiag/armadillo.hpp>

namespace xdiag::julia {

template <typename block_t> static void register_apply(jlcxx::Module &mod) {

  // Vector variants
  mod.method("fun_apply", [](OpSum const &ops, block_t const &bin, double *pin,
                             block_t const &bout, double *pout) {
    int64_t sin = bin.size();
    int64_t sout = bout.size();
    arma::vec vin(pin, sin, /*copy_aux_mem*/ false, /*strict*/ false);
    arma::vec vout(pout, sout, /*copy_aux_mem*/ false, /*strict*/ false);
    JULIA_XDIAG_CALL_VOID((apply(ops, bin, vin, bout, vout)));
  });

  mod.method("fun_applyC", [](OpSum const &ops, block_t const &bin,
                              complex *pin, block_t const &bout,
                              complex *pout) {
    int64_t sin = bin.size();
    int64_t sout = bout.size();
    arma::cx_vec vin(pin, sin, /*copy_aux_mem*/ false, /*strict*/ false);
    arma::cx_vec vout(pout, sout, /*copy_aux_mem*/ false, /*strict*/ false);
    JULIA_XDIAG_CALL_VOID((apply(ops, bin, vin, bout, vout)));
  });

  // Matrix variants
  mod.method("fun_apply", [](OpSum const &ops, block_t const &bin, double *pin,
                             block_t const &bout, double *pout, int64_t ncols) {
    int64_t sin = bin.size();
    int64_t sout = bout.size();
    arma::mat min(pin, sin, ncols, /*copy_aux_mem*/ false, /*strict*/ false);
    arma::mat mout(pout, sout, ncols, /*copy_aux_mem*/ false, /*strict*/ false);
    JULIA_XDIAG_CALL_VOID((apply(ops, bin, min, bout, mout)));
  });

  mod.method("fun_applyC", [](OpSum const &ops, block_t const &bin,
                              complex *pin, block_t const &bout, complex *pout,
                              int64_t ncols) {
    int64_t sin = bin.size();
    int64_t sout = bout.size();
    arma::cx_mat min(pin, sin, ncols, /*copy_aux_mem*/ false, /*strict*/ false);
    arma::cx_mat mout(pout, sout, ncols, /*copy_aux_mem*/ false,
                      /*strict*/ false);
    JULIA_XDIAG_CALL_VOID((apply(ops, bin, min, bout, mout)));
  });
}

void define_apply(jlcxx::Module &mod) {
  using namespace xdiag;
  register_apply<Spinhalf>(mod);
  register_apply<tJ>(mod);
  register_apply<Electron>(mod);
  register_apply<Boson>(mod);
  register_apply<Fermion>(mod);
}

} // namespace xdiag::julia
