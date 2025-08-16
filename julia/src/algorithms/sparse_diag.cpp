// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sparse_diag.hpp"

namespace xdiag::julia {

template <class block_t> static void define_eig0s(jlcxx::Module &mod) {
  mod.method("cxx_eigval0",
             [](OpSum const &ops, block_t const &block, double precision,
                int max_iterations, uint64_t seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigval0(ops, block, precision, max_iterations, seed));
             });

  mod.method("cxx_eig0", [](OpSum const &ops, block_t const &block,
                            double precision, int max_iterations,
                            uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eig0(ops, block, precision, max_iterations, seed));
  });
}

template <typename idx_t, typename coeff_t, class block_t>
static void define_eig0_sparse(jlcxx::Module &mod) {
  mod.method("cxx_eigval0",
             [](CSRMatrix<idx_t, coeff_t> const &ops, block_t const &block,
                double precision, int max_iterations, uint64_t seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigval0(ops, block, precision, max_iterations, seed));
             });

  mod.method("cxx_eig0", [](CSRMatrix<idx_t, coeff_t> const &ops,
                            block_t const &block, double precision,
                            int max_iterations, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eig0(ops, block, precision, max_iterations, seed));
  });
}

void define_eig0(jlcxx::Module &mod) {
  define_eig0s<Spinhalf>(mod);
  define_eig0s<tJ>(mod);
  define_eig0s<Electron>(mod);

  define_eig0_sparse<int64_t, double, Spinhalf>(mod);
  define_eig0_sparse<int64_t, double, tJ>(mod);
  define_eig0_sparse<int64_t, double, Electron>(mod);

  define_eig0_sparse<int64_t, complex, Spinhalf>(mod);
  define_eig0_sparse<int64_t, complex, tJ>(mod);
  define_eig0_sparse<int64_t, complex, Electron>(mod);

  define_eig0_sparse<int32_t, double, Spinhalf>(mod);
  define_eig0_sparse<int32_t, double, tJ>(mod);
  define_eig0_sparse<int32_t, double, Electron>(mod);

  define_eig0_sparse<int32_t, complex, Spinhalf>(mod);
  define_eig0_sparse<int32_t, complex, tJ>(mod);
  define_eig0_sparse<int32_t, complex, Electron>(mod);
}

} // namespace xdiag::julia
