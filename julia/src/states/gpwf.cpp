// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

namespace xdiag::julia {
void define_gpwf(jlcxx::Module &mod) {

  mod.add_type<GPWF>("cxx_GPWF")
      .constructor<>()
      .constructor<arma::mat const &, int64_t>()
      .constructor<arma::cx_mat const &, int64_t>()
      .method("nsites",
              [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("nup", [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.nup()) })
      .method("isreal",
              [](GPWF const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) });

  mod.method("to_string", [](GPWF const &s) { return to_string(s); });
}

} // namespace xdiag::julia
