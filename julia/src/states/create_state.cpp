// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

namespace xdiag::julia {

template <class block_t>
static void define_create_state_block(jlcxx::Module &mod) {
  mod.method("cxx_product_state",
             [](block_t const &block,
                std::vector<std::string> const &local_state, bool real) {
               JULIA_XDIAG_CALL_RETURN(product_state(block, local_state, real))
             });

  mod.method("cxx_random_state", [](block_t const &block, bool real,
                                    int64_t n_cols, int64_t seed,
                                    bool normalized) {
    JULIA_XDIAG_CALL_RETURN(
        random_state<block_t>(block, real, n_cols, seed, normalized));
  });

  mod.method("cxx_zero_state",
             [](block_t const &block, bool real, int64_t n_cols) {
               JULIA_XDIAG_CALL_RETURN(zero_state(block, real, n_cols));
             });
}

void define_create_state(jlcxx::Module &mod) {
  define_create_state_block<Spinhalf>(mod);
  define_create_state_block<tJ>(mod);
  define_create_state_block<Electron>(mod);

  mod.method("cxx_zero",
             [](State &state) { JULIA_XDIAG_CALL_VOID(zero(state)) });
}

} // namespace xdiag::julia
