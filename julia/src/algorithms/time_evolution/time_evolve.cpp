// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve.hpp"

namespace xdiag::julia {

template <typename op_t> static void define_time_evolve_op(jlcxx::Module &mod) {

  mod.method("cxx_time_evolve", [](op_t const &ops, State psi, double time,
                                   double precision, std::string algorithm) {
    JULIA_XDIAG_CALL_RETURN(time_evolve(ops, psi, time, precision, algorithm));
  });

  mod.method("cxx_time_evolve_inplace",
             [](op_t const &ops, State &psi, double time, double precision,
                std::string algorithm) {
               JULIA_XDIAG_CALL_RETURN(
                   time_evolve_inplace(ops, psi, time, precision, algorithm));
             });
}

void define_time_evolve(jlcxx::Module &mod) {
  define_time_evolve_op<OpSum>(mod);
  define_time_evolve_op<CSRMatrix<int64_t, double>>(mod);
  define_time_evolve_op<CSRMatrix<int64_t, complex>>(mod);
  define_time_evolve_op<CSRMatrix<int32_t, double>>(mod);
  define_time_evolve_op<CSRMatrix<int32_t, complex>>(mod);
} // namespace xdiag::julia
