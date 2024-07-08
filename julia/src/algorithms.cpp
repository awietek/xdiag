#include "algorithms.hpp"

namespace xdiag::julia {
    
void define_eigs_lanczos(jlcxx::Module &mod);
  
void define_eigval0(jlcxx::Module &mod) {
  // Computing the ground state energy
  mod.method("cxx_eigval0", [](BondList const &bonds, Spinhalf const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("cxx_eigval0", [](BondList const &bonds, tJ const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("cxx_eigval0", [](BondList const &bonds, Electron const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });
}

void define_eig0(jlcxx::Module &mod) {

  mod.method("cxx_eig0", [](BondList const &bonds, Spinhalf const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("cxx_eig0", [](BondList const &bonds, tJ const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("cxx_eig0", [](BondList const &bonds, Electron const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });
}

void define_exp_sym_v(jlcxx::Module &mod) {

  // methods for time evolution
  mod.method("cxx_exp_sym_v", [](BondList const &bonds, State state, double tau,
                                 bool normalize, double shift, double precision,
                                 int64_t max_iterations, double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(exp_sym_v(bonds, state, tau, normalize, shift,
                                      precision, max_iterations,
                                      deflation_tol));
  });

  mod.method("cxx_exp_sym_v", [](BondList const &bonds, State state,
                                 complex tau, bool normalize, double shift,
                                 double precision, int64_t max_iterations,
                                 double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(exp_sym_v(bonds, state, tau, normalize, shift,
                                      precision, max_iterations,
                                      deflation_tol));
  });
}

void define_time_evolve(jlcxx::Module &mod) {
  mod.method("cxx_time_evolve",
             [](BondList const &bonds, State state, double time,
                double precision, int64_t m, double anorm, int64_t nnorm) {
               JULIA_XDIAG_CALL_RETURN(
                   time_evolve(bonds, state, time, precision, m, anorm, nnorm));
             });
}

} // namespace xdiag::julia
