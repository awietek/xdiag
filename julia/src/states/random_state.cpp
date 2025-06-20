// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

namespace xdiag::julia {
void define_random_state(jlcxx::Module &mod) {
  mod.add_type<RandomState>("cxx_RandomState");

  mod.method("to_string", [](RandomState const &s) { return to_string(s); });

  mod.method("construct_RandomState", [](int64_t seed, bool normalized) {
    JULIA_XDIAG_CALL_RETURN(RandomState(seed, normalized))
  });
}

} // namespace xdiag::julia
