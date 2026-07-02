// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sparse_build.hpp"

#include <algorithm>
#include <numeric>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/blocks/dispatch_bases.hpp>
#include <xdiag/kernels/kernels.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

// Basis-level two-pass build, unified across CSR (transpose=false) and CSC
// (transpose=true). The transpose flag is forwarded to the shared csr kernels,
// which key the counting/fill by row (CSR) or column (CSC); everything else
// (offset prefix sum, ptr construction, per-group sort) is layout-agnostic.
template <typename idx_t, typename coeff_t, typename block_t, typename basis_t>
static void sparse_build_basis(OpSum const &ops, basis_t const &basis_in,
                               basis_t const &basis_out, idx_t ndim, idx_t i0,
                               bool transpose, arma::Col<idx_t> &ptr,
                               arma::Col<idx_t> &idx, arma::Col<coeff_t> &data) {
  // Pass 1: count nonzeros per group (row for CSR, column for CSC).
  auto counts = kernels::csr_matrix_nnz<block_t, coeff_t>(ops, basis_in,
                                                          basis_out, transpose);
  int64_t nnz = std::accumulate(counts.begin(), counts.end(), (int64_t)0);

  // Build ptr (exclusive prefix sum + i0 shift) and the mutable offset array
  // that csr_matrix_fill uses to claim write slots.
  ptr.resize(ndim + 1);
  idx.resize(nnz);
  data.resize(nnz);

  std::vector<int64_t> offset(ndim, 0);
  if (ndim > 0) {
    std::partial_sum(counts.begin(), counts.end() - 1, offset.begin() + 1);
    for (idx_t k = 0; k < ndim; ++k) {
      ptr[k] = (idx_t)(offset[k] + i0);
    }
    ptr[ndim] = (idx_t)(nnz + i0);
  }

  // Pass 2: fill idx and data using atomic slot assignment.
  kernels::csr_matrix_fill<block_t, coeff_t>(ops, basis_in, basis_out, offset,
                                             idx.memptr(), data.memptr(), i0,
                                             transpose);

  // Sort each group's entries by stored index (required for CSR/CSC validity).
  int64_t max_elems =
      counts.empty() ? 0 : *std::max_element(counts.begin(), counts.end());
#ifdef _OPENMP
#pragma omp parallel
  {
    std::vector<int64_t> indices(max_elems);
    std::vector<idx_t> idxtmp(max_elems);
    std::vector<coeff_t> datatmp(max_elems);
#pragma omp for
#else
  {
    std::vector<int64_t> indices(max_elems);
    std::vector<idx_t> idxtmp(max_elems);
    std::vector<coeff_t> datatmp(max_elems);
#endif
    for (idx_t g = 0; g < ndim; ++g) {
      int64_t start = (int64_t)(ptr[g] - i0);
      int64_t end = (int64_t)(ptr[g + 1] - i0);
      int64_t nelems = end - start;
      std::copy(idx.memptr() + start, idx.memptr() + end, idxtmp.begin());
      std::copy(data.memptr() + start, data.memptr() + end, datatmp.begin());
      std::iota(indices.begin(), indices.begin() + nelems, (int64_t)0);
      std::sort(indices.begin(), indices.begin() + nelems,
                [&](int64_t i, int64_t j) { return idxtmp[i] < idxtmp[j]; });
      for (int64_t i = 0; i < nelems; ++i) {
        idx[start + i] = idxtmp[indices[i]];
        data[start + i] = datatmp[indices[i]];
      }
    }
  }
}

template <typename idx_t, typename coeff_t, typename block_t>
void build_csr_arrays(OpSum const &ops, block_t const &block_in,
                      block_t const &block_out, idx_t ndim, idx_t i0,
                      bool transpose, arma::Col<idx_t> &ptr,
                      arma::Col<idx_t> &idx, arma::Col<coeff_t> &data) try {
  kernels::dispatch_basis(
      block_in, block_out, [&](auto const &basis_in, auto const &basis_out) {
        sparse_build_basis<idx_t, coeff_t, block_t>(
            ops, basis_in, basis_out, ndim, i0, transpose, ptr, idx, data);
      });
}
XDIAG_CATCH

// Explicit instantiations for every (idx_t, coeff_t, non-distributed block).
// This is the ONLY place kernels::dispatch_basis is instantiated for the sparse
// build, so csr_matrix.cpp / csc_matrix.cpp share one copy of the visitor.
#define XDIAG_INST_BUILD_IC(BLOCK, IDX, COEFF)                                 \
  template void build_csr_arrays<IDX, COEFF, BLOCK>(                           \
      OpSum const &, BLOCK const &, BLOCK const &, IDX, IDX, bool,             \
      arma::Col<IDX> &, arma::Col<IDX> &, arma::Col<COEFF> &);

#define XDIAG_INST_BUILD(BLOCK)                                                \
  XDIAG_INST_BUILD_IC(BLOCK, int32_t, double)                                  \
  XDIAG_INST_BUILD_IC(BLOCK, int32_t, complex)                                 \
  XDIAG_INST_BUILD_IC(BLOCK, int64_t, double)                                  \
  XDIAG_INST_BUILD_IC(BLOCK, int64_t, complex)

XDIAG_INST_BUILD(Spinhalf)
XDIAG_INST_BUILD(Boson)
XDIAG_INST_BUILD(Fermion)
XDIAG_INST_BUILD(Electron)
XDIAG_INST_BUILD(tJ)

#undef XDIAG_INST_BUILD
#undef XDIAG_INST_BUILD_IC

} // namespace xdiag
