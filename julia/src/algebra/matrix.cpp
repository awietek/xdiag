#include "matrix.hpp"

#include <xdiag/all.hpp>

namespace xdiag::julia {

template <typename coeff_t, class block_t>
static void define_matrices(jlcxx::Module &mod) {
  std::string name = isreal<coeff_t>() ? "cxx_matrix" : "cxx_matrixC";

  mod.method(name, [](coeff_t *mat, Op const &op, block_t const &block) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, op, block));
  });
  mod.method(name, [](coeff_t *mat, OpSum const &ops, block_t const &block) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block));
  });

  mod.method(name, [](coeff_t *mat, Op const &op, block_t const &block_in,
                      block_t const &block_out) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, op, block_in, block_out));
  });
  mod.method(name, [](coeff_t *mat, OpSum const &ops, block_t const &block_in,
                      block_t const &block_out) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out));
  });

}

void define_matrix(jlcxx::Module &mod) {
  define_matrices<double, Spinhalf>(mod);
  define_matrices<double, tJ>(mod);
  define_matrices<double, Electron>(mod);

  define_matrices<complex, Spinhalf>(mod);
  define_matrices<complex, tJ>(mod);
  define_matrices<complex, Electron>(mod);
}

} // namespace xdiag::julia
