// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "opsum.hpp"

namespace xdiag::julia {

void define_opsum(jlcxx::Module &mod) {

  mod.add_type<OpSum>("cxx_OpSum")
      .method("plain",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(ops.plain()); })
      .method("isreal",
              [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(isreal(ops)); })
      .method("+", [](OpSum const &ops,
                      Op const &op) { JULIA_XDIAG_CALL_RETURN(ops + op); })
      .method("+",
              [](OpSum const &ops, OpSum const &ops2) {
                JULIA_XDIAG_CALL_RETURN(ops + ops2)
              })
      .method("-", [](OpSum const &ops,
                      Op const &op) { JULIA_XDIAG_CALL_RETURN(ops - op); })
      .method("-",
              [](OpSum const &ops, OpSum const &ops2) {
                JULIA_XDIAG_CALL_RETURN(ops - ops2);
              })
      .method("*", [](int64_t cpl,
                      OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](double cpl,
                      OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](complex cpl,
                      OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](OpSum const &ops,
                      int64_t cpl) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](OpSum const &ops,
                      double cpl) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("*", [](OpSum const &ops,
                      complex cpl) { JULIA_XDIAG_CALL_RETURN(cpl * ops); })
      .method("/", [](OpSum const &ops,
                      int64_t cpl) { JULIA_XDIAG_CALL_RETURN(ops / cpl); })
      .method("/", [](OpSum const &ops,
                      double cpl) { JULIA_XDIAG_CALL_RETURN(ops / cpl); })
      .method("/", [](OpSum const &ops,
                      complex cpl) { JULIA_XDIAG_CALL_RETURN(ops / cpl); })
      .method("setindex!",
              [](OpSum &ops, std::string name, int64_t value) {
                JULIA_XDIAG_CALL_VOID(ops[name] = (double)value);
              })
      .method("setindex!",
              [](OpSum &ops, std::string name, double value) {
                JULIA_XDIAG_CALL_VOID(ops[name] = value);
              })
      .method("setindex!",
              [](OpSum &ops, std::string name, complex value) {
                JULIA_XDIAG_CALL_VOID(ops[name] = value);
              })
      .method("constants", [](OpSum const &ops) {
        JULIA_XDIAG_CALL_RETURN(ops.constants());
      });

  mod.method("==", [](OpSum const &a, OpSum const &b) {
    JULIA_XDIAG_CALL_RETURN(a == b);
  });
  mod.method("to_string",
             [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(to_string(ops)); });
  mod.method("isapprox", [](OpSum const &ops1, OpSum const &ops2, double rtol,
                            double atol) {
    JULIA_XDIAG_CALL_RETURN(isapprox(ops1, ops2, rtol, atol));
  });

  mod.method("construct_OpSum", []() { JULIA_XDIAG_CALL_RETURN(OpSum()); });
  mod.method("construct_OpSum",
             [](Op const &op) { JULIA_XDIAG_CALL_RETURN(OpSum(op)); });
  mod.method("construct_OpSum", [](std::string cpl, Op const &op) {
    JULIA_XDIAG_CALL_RETURN(OpSum(cpl, op));
  });
  mod.method("construct_OpSum", [](int64_t cpl, Op const &op) {
    JULIA_XDIAG_CALL_RETURN(OpSum(cpl, op));
  });
  mod.method("construct_OpSum", [](double cpl, Op const &op) {
    JULIA_XDIAG_CALL_RETURN(OpSum(cpl, op));
  });
  mod.method("construct_OpSum", [](complex cpl, Op const &op) {
    JULIA_XDIAG_CALL_RETURN(OpSum(cpl, op));
  });
}
} // namespace xdiag::julia
