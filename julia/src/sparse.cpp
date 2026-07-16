// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// The output block is inferred here via block(ops, in), so Julia never has to
// round-trip the Block variant. Phase-1 (nnz) and phase-2 (fill) both re-infer
// it (deterministic). Extra parens around the templated calls keep the comma
// in <idx,coeff> from splitting the JULIA_XDIAG_CALL_* macro arguments.

#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {

template <typename block_t> static void register_sparse(jlcxx::Module &mod) {
  // number of rows = size of the (inferred) output block
  mod.method("fun_sparse_dim_out", [](OpSum const &ops, block_t const &in) {
    JULIA_XDIAG_CALL_RETURN((int64_t)size(block(ops, in)));
  });

  // --- COO -----------------------------------------------------------------
  // nnz is structural but the kernel enforces coeff_t == realness(ops, block).
  mod.method("fun_coo_nnz", [](OpSum const &ops, block_t const &in) {
    if (isreal(ops) && isreal(in)) {
      JULIA_XDIAG_CALL_RETURN(coo_matrix_nnz<double>(ops, in, block(ops, in)));
    } else {
      JULIA_XDIAG_CALL_RETURN(coo_matrix_nnz<complex>(ops, in, block(ops, in)));
    }
  });
  mod.method("fun_coo_fill",
             [](OpSum const &ops, block_t const &in, int64_t nnz, int64_t *row,
                int64_t *col, double *data, int64_t i0) {
               JULIA_XDIAG_CALL_VOID((coo_matrix_fill<int64_t, double>(
                   ops, in, block(ops, in), nnz, row, col, data, i0)));
             });
  mod.method("fun_coo_fill",
             [](OpSum const &ops, block_t const &in, int64_t nnz, int32_t *row,
                int32_t *col, double *data, int32_t i0) {
               JULIA_XDIAG_CALL_VOID((coo_matrix_fill<int32_t, double>(
                   ops, in, block(ops, in), nnz, row, col, data, i0)));
             });
  mod.method("fun_coo_fill",
             [](OpSum const &ops, block_t const &in, int64_t nnz, int64_t *row,
                int64_t *col, complex *data, int64_t i0) {
               JULIA_XDIAG_CALL_VOID((coo_matrix_fill<int64_t, complex>(
                   ops, in, block(ops, in), nnz, row, col, data, i0)));
             });
  mod.method("fun_coo_fill",
             [](OpSum const &ops, block_t const &in, int64_t nnz, int32_t *row,
                int32_t *col, complex *data, int32_t i0) {
               JULIA_XDIAG_CALL_VOID((coo_matrix_fill<int32_t, complex>(
                   ops, in, block(ops, in), nnz, row, col, data, i0)));
             });

  // --- CSR / CSC (transpose=true builds CSC arrays directly) ---------------
  mod.method(
      "fun_csr_nnz", [](OpSum const &ops, block_t const &in, bool transpose) {
        if (isreal(ops) && isreal(in)) {
          JULIA_XDIAG_CALL_RETURN(
              csr_matrix_nnz<double>(ops, in, block(ops, in), transpose));
        } else {
          JULIA_XDIAG_CALL_RETURN(
              csr_matrix_nnz<complex>(ops, in, block(ops, in), transpose));
        }
      });
  mod.method("fun_csr_fill", [](OpSum const &ops, block_t const &in,
                                std::vector<int64_t> const &counts,
                                int64_t *ptr, int64_t *idx, double *data,
                                int64_t i0, bool transpose) {
    JULIA_XDIAG_CALL_VOID((csr_matrix_fill<int64_t, double>(
        ops, in, block(ops, in), counts, ptr, idx, data, i0, transpose)));
  });
  mod.method("fun_csr_fill", [](OpSum const &ops, block_t const &in,
                                std::vector<int64_t> const &counts,
                                int32_t *ptr, int32_t *idx, double *data,
                                int32_t i0, bool transpose) {
    JULIA_XDIAG_CALL_VOID((csr_matrix_fill<int32_t, double>(
        ops, in, block(ops, in), counts, ptr, idx, data, i0, transpose)));
  });
  mod.method("fun_csr_fill", [](OpSum const &ops, block_t const &in,
                                std::vector<int64_t> const &counts,
                                int64_t *ptr, int64_t *idx, complex *data,
                                int64_t i0, bool transpose) {
    JULIA_XDIAG_CALL_VOID((csr_matrix_fill<int64_t, complex>(
        ops, in, block(ops, in), counts, ptr, idx, data, i0, transpose)));
  });
  mod.method("fun_csr_fill", [](OpSum const &ops, block_t const &in,
                                std::vector<int64_t> const &counts,
                                int32_t *ptr, int32_t *idx, complex *data,
                                int32_t i0, bool transpose) {
    JULIA_XDIAG_CALL_VOID((csr_matrix_fill<int32_t, complex>(
        ops, in, block(ops, in), counts, ptr, idx, data, i0, transpose)));
  });
}

void define_sparse(jlcxx::Module &mod) {
  using namespace xdiag;
  register_sparse<Spinhalf>(mod);
  register_sparse<tJ>(mod);
  register_sparse<Electron>(mod);
  register_sparse<Boson>(mod);
  register_sparse<Fermion>(mod);
}

} // namespace xdiag::julia
