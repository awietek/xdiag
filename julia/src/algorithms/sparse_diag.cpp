#include "sparse_diag.hpp"

namespace xdiag::julia {

template <class block_t> static void define_eig0s(jlcxx::Module &mod) {
  mod.method("eigval0",
             [](OpSum const &ops, block_t const &block, double precision,
                int max_iterations, uint64_t seed) {
               JULIA_XDIAG_CALL_RETURN(
                   eigval0(ops, block, precision, max_iterations, seed));
             });

  mod.method("eig0", [](OpSum const &ops, block_t const &block,
                        double precision, int max_iterations, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(eig0(ops, block, precision, max_iterations, seed));
  });
}

void define_eig0(jlcxx::Module &mod) {
  define_eig0s<Spinhalf>(mod);
  define_eig0s<tJ>(mod);
  define_eig0s<Electron>(mod);
}

} // namespace xdiag::julia
