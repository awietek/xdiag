#include "sparse_diag.hpp"

namespace xdiag::julia {

void define_eigval0(jlcxx::Module &mod) {
  // Computing the ground state energy
  mod.method("eigval0", [](OpSum const &ops, Spinhalf const &block,
                           double precision, int max_iterations,
                           bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(ops, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eigval0", [](OpSum const &ops, tJ const &block, double precision,
                           int max_iterations, bool force_complex,
                           uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(ops, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eigval0", [](OpSum const &ops, Electron const &block,
                           double precision, int max_iterations,
                           bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(ops, block, precision, max_iterations, force_complex, seed));
  });
}

void define_eig0(jlcxx::Module &mod) {

  mod.method("eig0", [](OpSum const &ops, Spinhalf const &block,
                        double precision, int max_iterations,
                        bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(ops, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0", [](OpSum const &ops, tJ const &block, double precision,
                        int max_iterations, bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(ops, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0", [](OpSum const &ops, Electron const &block,
                        double precision, int max_iterations,
                        bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(ops, block, precision, max_iterations, force_complex, seed));
  });
}

} // namespace xdiag::julia
