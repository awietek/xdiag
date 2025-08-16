// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "time_evolve_expokit.hpp"

namespace xdiag::julia {

template <typename op_t>
static void define_time_evolve_expokit_op(jlcxx::Module &mod) {
  mod.method("cxx_time_evolve_expokit",
             [](op_t const &H, State psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit(
                   H, psi, time, precision, m, anorm, nnorm));
             });

  mod.method("cxx_time_evolve_expokit_inplace",
             [](op_t const &H, State &psi, double time, double precision,
                int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(time_evolve_expokit_inplace(
                   H, psi, time, precision, m, anorm, nnorm));
             });
}

void define_time_evolve_expokit(jlcxx::Module &mod) {

  using res_t = TimeEvolveExpokitResult;
  mod.add_type<res_t>("cxx_TimeEvolveExpokitResult")
      .method("error",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.error) })
      .method("hump",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.hump) })
      .method("state",
              [](res_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.state) });

  using res_inplace_t = TimeEvolveExpokitInplaceResult;
  mod.add_type<res_inplace_t>("cxx_TimeEvolveExpokitInplaceResult")
      .method(
          "error",
          [](res_inplace_t const &r) { JULIA_XDIAG_CALL_RETURN_MOVE(r.error) })
      .method("hump", [](res_inplace_t const &r) {
        JULIA_XDIAG_CALL_RETURN_MOVE(r.hump)
      });
  define_time_evolve_expokit_op<OpSum>(mod);
  define_time_evolve_expokit_op<CSRMatrix<int64_t, double>>(mod);
  define_time_evolve_expokit_op<CSRMatrix<int64_t, complex>>(mod);
  define_time_evolve_expokit_op<CSRMatrix<int32_t, double>>(mod);
  define_time_evolve_expokit_op<CSRMatrix<int32_t, complex>>(mod);
}
} // namespace xdiag::julia
