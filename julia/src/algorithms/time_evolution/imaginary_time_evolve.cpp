// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "imaginary_time_evolve.hpp"

namespace xdiag::julia {

template <typename op_t>
static void define_imaginary_time_evolve_op(jlcxx::Module &mod) {

  mod.method("cxx_imaginary_time_evolve",
             [](op_t const &ops, State psi, double time, double precision,
                double shift) {
               JULIA_XDIAG_CALL_RETURN(
                   imaginary_time_evolve(ops, psi, time, precision, shift));
             });

  mod.method("cxx_imaginary_time_evolve_inplace",
             [](op_t const &ops, State &psi, double time, double precision,
                double shift) {
               JULIA_XDIAG_CALL_RETURN(imaginary_time_evolve_inplace(
                   ops, psi, time, precision, shift));
             });
}

void define_imaginary_time_evolve(jlcxx::Module &mod) {
  define_imaginary_time_evolve_op<OpSum>(mod);
  define_imaginary_time_evolve_op<CSRMatrix<int64_t, double>>(mod);
  define_imaginary_time_evolve_op<CSRMatrix<int64_t, complex>>(mod);
  define_imaginary_time_evolve_op<CSRMatrix<int32_t, double>>(mod);
  define_imaginary_time_evolve_op<CSRMatrix<int32_t, complex>>(mod);
}
} // namespace xdiag::julia
