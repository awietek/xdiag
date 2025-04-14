// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply.hpp"
#include <xdiag/all.hpp>

namespace xdiag::julia {

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
}

} // namespace xdiag::julia
