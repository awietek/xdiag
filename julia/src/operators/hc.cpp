// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "hc.hpp"

namespace xdiag::julia {

void define_hc(jlcxx::Module &mod) {

  mod.method("cxx_hc", [](Op const &op) { JULIA_XDIAG_CALL_RETURN(hc(op)); });
  mod.method("cxx_hc",
             [](OpSum const &ops) { JULIA_XDIAG_CALL_RETURN(hc(ops)); });
  mod.method("cxx_ishermitian", [](OpSum const &ops) {
    JULIA_XDIAG_CALL_RETURN(ishermitian(ops));
  });
}

} // namespace xdiag::julia
