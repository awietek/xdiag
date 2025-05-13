// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {

void define_algebra(jlcxx::Module &mod) {

  mod.method("cxx_norm",
             [](State const &v) { JULIA_XDIAG_CALL_RETURN(norm(v)); });
  mod.method("cxx_norm1",
             [](State const &v) { JULIA_XDIAG_CALL_RETURN(norm1(v)); });
  mod.method("cxx_norminf",
             [](State const &v) { JULIA_XDIAG_CALL_RETURN(norminf(v)); });

  mod.method("cxx_dot", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(dot(v, w));
  });

  mod.method("cxx_dotC", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(dotC(v, w));
  });

  mod.method("cxx_inner", [](OpSum const &ops, State const &v) {
    JULIA_XDIAG_CALL_RETURN(inner(ops, v));
  });

  mod.method("cxx_inner", [](Op const &op, State const &v) {
    JULIA_XDIAG_CALL_RETURN(inner(op, v));
  });

  mod.method("cxx_innerC", [](OpSum const &ops, State const &v) {
    JULIA_XDIAG_CALL_RETURN(innerC(ops, v));
  });

  mod.method("cxx_innerC", [](Op const &op, State const &v) {
    JULIA_XDIAG_CALL_RETURN(innerC(op, v));
  });

  mod.method("cxx_add", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(v + w);
  });

  mod.method("cxx_sub", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(v - w);
  });


  mod.method("cxx_scalar_mult", [](double alpha, State const &v) {
    JULIA_XDIAG_CALL_RETURN(alpha * v);
  });

  mod.method("cxx_scalar_mult", [](complex alpha, State const &v) {
    JULIA_XDIAG_CALL_RETURN(alpha * v);
  });

  mod.method("cxx_scalar_div", [](State const &v, double alpha) {
    JULIA_XDIAG_CALL_RETURN(v / alpha);
  });

  mod.method("cxx_scalar_div", [](State const &v, complex alpha) {
    JULIA_XDIAG_CALL_RETURN(v / alpha);
  });
}

} // namespace xdiag::julia
