// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "op.hpp"

namespace xdiag::julia {

void define_op(jlcxx::Module &mod) {

  mod.add_type<Op>("cxx_Op")
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

  mod.method("construct_Op", []() { JULIA_XDIAG_CALL_RETURN(Op()); });
  mod.method("construct_Op",
             [](std::string type) { JULIA_XDIAG_CALL_RETURN(Op(type)); });
  mod.method("construct_Op", [](std::string type, int64_t site) {
    JULIA_XDIAG_CALL_RETURN(Op(type, site));
  });
  mod.method("construct_Op",
             [](std::string type, std::vector<int64_t> const &sites) {
               JULIA_XDIAG_CALL_RETURN(Op(type, sites));
             });
  mod.method("construct_Op",
             [](std::string type, int64_t site, arma::mat const &mat) {
               JULIA_XDIAG_CALL_RETURN(Op(type, site, mat));
             });
  mod.method("construct_Op",
             [](std::string type, std::vector<int64_t> const &sites,
                arma::mat const &mat) {
               JULIA_XDIAG_CALL_RETURN(Op(type, sites, mat));
             });
  mod.method("construct_Op",
             [](std::string type, int64_t site, arma::cx_mat const &mat) {
               JULIA_XDIAG_CALL_RETURN(Op(type, site, mat));
             });
  mod.method("construct_Op",
             [](std::string type, std::vector<int64_t> const &sites,
                arma::cx_mat const &mat) {
               JULIA_XDIAG_CALL_RETURN(Op(type, sites, mat));
             });
}

} // namespace xdiag::julia
