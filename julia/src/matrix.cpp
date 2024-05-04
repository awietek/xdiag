#include "matrix.hpp"

namespace xdiag::julia {

void define_matrix_cxx(jlcxx::Module &mod) {

  mod.method("matrix_cxx",
             [](double *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
             });

  mod.method("matrixC_cxx",
             [](complex *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrixC(reinterpret_cast<complex *>(mat),
                                             bonds, block_in, block_out));
             });

  mod.method("matrix_cxx", [](double *mat, BondList const &bonds,
                              tJ const &block_in, tJ const &block_out) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
  });

  mod.method("matrixC_cxx", [](complex *mat, BondList const &bonds,
                               tJ const &block_in, tJ const &block_out) {
    JULIA_XDIAG_CALL_VOID(
        matrixC(reinterpret_cast<complex *>(mat), bonds, block_in, block_out));
  });

  mod.method("matrix_cxx",
             [](double *mat, BondList const &bonds, Electron const &block_in,
                Electron const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
             });

  mod.method("matrixC_cxx",
             [](complex *mat, BondList const &bonds, Electron const &block_in,
                Electron const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrixC(reinterpret_cast<complex *>(mat),
                                             bonds, block_in, block_out));
             });
}

} // namespace xdiag::julia
