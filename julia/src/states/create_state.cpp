#include "state.hpp"

namespace xdiag::julia {
void define_create_state(jlcxx::Module &mod) {

  mod.method("cxx_product",
             [](Spinhalf const &block,
                std::vector<std::string> const &local_state, bool real) {
               JULIA_XDIAG_CALL_RETURN(product(block, local_state, real))
             });
  mod.method("cxx_product",
             [](tJ const &block, std::vector<std::string> const &local_state,
                bool real) {
               JULIA_XDIAG_CALL_RETURN(product(block, local_state, real))
             });
  mod.method("cxx_product",
             [](Electron const &block,
                std::vector<std::string> const &local_state, bool real) {
               JULIA_XDIAG_CALL_RETURN(product(block, local_state, real))
             });

  mod.method("cxx_rand", [](Spinhalf const &block, bool real, int64_t seed,
                            bool normalized) {
    JULIA_XDIAG_CALL_RETURN(rand(block, real, seed, normalized))
  });
  mod.method("cxx_rand",
             [](tJ const &block, bool real, int64_t seed, bool normalized) {
               JULIA_XDIAG_CALL_RETURN(rand(block, real, seed, normalized))
             });
  mod.method("cxx_rand", [](Electron const &block, bool real, int64_t seed,
                            bool normalized) {
    JULIA_XDIAG_CALL_RETURN(rand(block, real, seed, normalized))
  });

  mod.method("cxx_zeros", [](Spinhalf const &block, bool real, int64_t n_cols) {
    JULIA_XDIAG_CALL_RETURN(zeros(block, real, n_cols))
  });
  mod.method("cxx_zeros", [](tJ const &block, bool real, int64_t n_cols) {
    JULIA_XDIAG_CALL_RETURN(zeros(block, real, n_cols))
  });
  mod.method("cxx_zeros", [](Electron const &block, bool real, int64_t n_cols) {
    JULIA_XDIAG_CALL_RETURN(zeros(block, real, n_cols))
  });

  mod.method("cxx_zero",
             [](State &state) { JULIA_XDIAG_CALL_VOID(zero(state)) });
}

} // namespace xdiag::julia
