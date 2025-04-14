// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op.hpp"

namespace xdiag::julia {

void define_op(jlcxx::Module &mod) {

  mod.add_type<Op>("cxx_Op")
      .constructor<>()
      .constructor<std::string>()

      .constructor<std::string, int64_t>()
      .constructor<std::string, std::vector<int64_t> const &>()

      .constructor<std::string, int64_t, arma::mat const &>()
      .constructor<std::string, std::vector<int64_t> const &,
                   arma::mat const &>()

      .constructor<std::string, int64_t, arma::cx_mat const &>()
      .constructor<std::string, std::vector<int64_t> const &,
                   arma::cx_mat const &>()

      .method("type", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.type()) })
      .method("size", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.size()) })
      .method("getindex",
              [](Op const &op, int64_t idx) {
                JULIA_XDIAG_CALL_RETURN(op[idx - 1])
              })
      .method("sites", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(op.sites()) })
      .method("isreal",
              [](Op const &op) { JULIA_XDIAG_CALL_RETURN(isreal(op)) });

  mod.method("==", [](Op const &a, Op const &b) { return a == b; });
  mod.method("to_string", [](Op const &op) { return to_string(op); });
  mod.method("isapprox",
             [](Op const &op1, Op const &op2, double rtol, double atol) {
               return isapprox(op1, op2, rtol, atol);
             });
}

} // namespace xdiag::julia
