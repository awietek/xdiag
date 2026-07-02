// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
//
// Hand-written specials for the Julia wrapper: the pointer-fill routines that
// let Julia allocate the output storage and have C++ fill it in place (no
// copy). These cannot be produced by the mechanical emitter (they take raw
// Julia-owned pointers and are templated over idx/coeff/block), so they are
// maintained here and registered from the generated module via define_specials.
//
//   dense:  matrix<coeff_t>(ops, block_in, block_out, coeff_t *mat)
//   sparse: csr_matrix_nnz  -> per-row counts (Julia sums + allocates)
//           csr_matrix_fill -> fills Julia-owned rowptr/col/data
//
// Block-variant parameters are expanded over the concrete serial blocks (CxxWrap
// cannot marshal std::variant); the concrete block converts implicitly to Block.

#include <julia/src/xdiagjl.hpp>

#include <xdiag/kernels/matrix.hpp>
#include <xdiag/kernels/sparse/coo_matrix.hpp>
#include <xdiag/kernels/sparse/csr_matrix.hpp>

namespace xdiag::julia {

// --- dense pointer-fill --------------------------------------------------
template <typename coeff_t, typename block_t>
static void define_matrix_fill(jlcxx::Module &mod) {
  std::string name =
      std::is_same_v<coeff_t, double> ? "cxx_matrix_fill" : "cxx_matrixC_fill";
  mod.method(name, [](OpSum const &ops, block_t const &block_in,
                      block_t const &block_out, coeff_t *mat) {
    JULIA_XDIAG_CALL_VOID(matrix(ops, block_in, block_out, mat));
  });
}

// --- sparse CSR two-phase -----------------------------------------------
template <typename coeff_t, typename block_t>
static void define_csr_nnz(jlcxx::Module &mod) {
  std::string name = std::is_same_v<coeff_t, double> ? "cxx_csr_matrix_nnz"
                                                     : "cxx_csr_matrixC_nnz";
  mod.method(name, [](OpSum const &ops, block_t const &block_in,
                      block_t const &block_out) {
    JULIA_XDIAG_CALL_RETURN(csr_matrix_nnz<coeff_t>(ops, block_in, block_out));
  });
}

template <typename idx_t, typename coeff_t, typename block_t>
static void define_csr_fill(jlcxx::Module &mod, std::string const &name) {
  mod.method(name, [](OpSum const &ops, block_t const &block_in,
                      block_t const &block_out,
                      std::vector<int64_t> const &n_elements_in_row,
                      idx_t *rowptr, idx_t *col, coeff_t *data, idx_t i0) {
    // idx_t/coeff_t are deduced from the pointer args (explicit template args
    // would put a comma inside the JULIA_XDIAG_CALL_VOID macro).
    JULIA_XDIAG_CALL_VOID(csr_matrix_fill(ops, block_in, block_out,
                                          n_elements_in_row, rowptr, col, data,
                                          i0));
  });
}

// --- sparse COO two-phase -----------------------------------------------
template <typename coeff_t, typename block_t>
static void define_coo_nnz(jlcxx::Module &mod) {
  std::string name = std::is_same_v<coeff_t, double> ? "cxx_coo_matrix_nnz"
                                                     : "cxx_coo_matrixC_nnz";
  mod.method(name, [](OpSum const &ops, block_t const &block_in,
                      block_t const &block_out) {
    JULIA_XDIAG_CALL_RETURN(coo_matrix_nnz<coeff_t>(ops, block_in, block_out));
  });
}

template <typename idx_t, typename coeff_t, typename block_t>
static void define_coo_fill(jlcxx::Module &mod, std::string const &name) {
  mod.method(name, [](OpSum const &ops, block_t const &block_in,
                      block_t const &block_out, int64_t nnz_capacity,
                      idx_t *row, idx_t *col, coeff_t *data, idx_t i0) {
    JULIA_XDIAG_CALL_VOID(coo_matrix_fill(ops, block_in, block_out,
                                          nnz_capacity, row, col, data, i0));
  });
}

template <typename block_t> static void define_block_specials(jlcxx::Module &mod) {
  define_matrix_fill<double, block_t>(mod);
  define_matrix_fill<complex, block_t>(mod);

  define_csr_nnz<double, block_t>(mod);
  define_csr_nnz<complex, block_t>(mod);
  define_csr_fill<int64_t, double, block_t>(mod, "cxx_csr_matrix_fill");
  define_csr_fill<int64_t, complex, block_t>(mod, "cxx_csr_matrixC_fill");
  define_csr_fill<int32_t, double, block_t>(mod, "cxx_csr_matrix_32_fill");
  define_csr_fill<int32_t, complex, block_t>(mod, "cxx_csr_matrixC_32_fill");

  define_coo_nnz<double, block_t>(mod);
  define_coo_nnz<complex, block_t>(mod);
  define_coo_fill<int64_t, double, block_t>(mod, "cxx_coo_matrix_fill");
  define_coo_fill<int64_t, complex, block_t>(mod, "cxx_coo_matrixC_fill");
  define_coo_fill<int32_t, double, block_t>(mod, "cxx_coo_matrix_32_fill");
  define_coo_fill<int32_t, complex, block_t>(mod, "cxx_coo_matrixC_32_fill");
}

void define_specials(jlcxx::Module &mod) {
  define_block_specials<Spinhalf>(mod);
  define_block_specials<tJ>(mod);
  define_block_specials<Electron>(mod);
  define_block_specials<Boson>(mod);
  define_block_specials<Fermion>(mod);
}

} // namespace xdiag::julia
