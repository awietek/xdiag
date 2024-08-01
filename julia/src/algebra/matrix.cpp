#include "matrix.hpp"

namespace xdiag::julia {

void define_matrix(jlcxx::Module &mod) {

  // Spinhalf
  mod.method("matrix", [](double *mat, OpSum const &ops,
                          Spinhalf const &block_in, Spinhalf const &block_out,
                          double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });

  mod.method("matrixC", [](complex *mat, OpSum const &ops,
                           Spinhalf const &block_in, Spinhalf const &block_out,
                           double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });

  // tJ
  mod.method("matrix", [](double *mat, OpSum const &ops, tJ const &block_in,
                          tJ const &block_out, double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });

  mod.method("matrix", [](complex *mat, OpSum const &ops, tJ const &block_in,
                          tJ const &block_out, double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });

  // Electron
  mod.method("matrix", [](double *mat, OpSum const &ops,
                          Electron const &block_in, Electron const &block_out,
                          double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });

  mod.method("matrix", [](complex *mat, OpSum const &ops,
                          Electron const &block_in, Electron const &block_out,
                          double precision) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, ops, block_in, block_out, precision));
  });
}

} // namespace xdiag::julia
