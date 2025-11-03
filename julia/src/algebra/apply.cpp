// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename block_t, typename mat_t>
static void define_apply_block(jlcxx::Module &mod) {

  mod.method("cxx_apply", [](OpSum const &ops, block_t const &block_in,
                             mat_t const &mat_in, block_t block_out,
                             mat_t &mat_out) {
    JULIA_XDIAG_CALL_VOID(apply(ops, block_in, mat_in, block_out, mat_out));
  });
}

void define_apply(jlcxx::Module &mod) {

  mod.method("cxx_apply", [](Op const &op, State const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(op, v));
  });

  mod.method("cxx_apply", [](OpSum const &ops, State const &v) {
    JULIA_XDIAG_CALL_RETURN(apply(ops, v));
  });

  mod.method("cxx_apply", [](Op const &op, State const &v, State &w) {
    JULIA_XDIAG_CALL_VOID(apply(op, v, w));
  });

  mod.method("cxx_apply", [](OpSum const &ops, State const &v, State &w) {
    JULIA_XDIAG_CALL_VOID(apply(ops, v, w));
  });

  define_apply_block<Spinhalf, arma::vec>(mod);
  define_apply_block<Spinhalf, arma::cx_vec>(mod);
  define_apply_block<Spinhalf, arma::mat>(mod);
  define_apply_block<Spinhalf, arma::cx_mat>(mod);

  define_apply_block<tJ, arma::vec>(mod);
  define_apply_block<tJ, arma::cx_vec>(mod);
  define_apply_block<tJ, arma::mat>(mod);
  define_apply_block<tJ, arma::cx_mat>(mod);

  define_apply_block<Electron, arma::vec>(mod);
  define_apply_block<Electron, arma::cx_vec>(mod);
  define_apply_block<Electron, arma::mat>(mod);
  define_apply_block<Electron, arma::cx_mat>(mod);
}

} // namespace xdiag::julia
