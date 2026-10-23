// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "sparse_build.hpp"

#include <algorithm>
#include <numeric>

#include <xdiag/kernels/elements.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/utils/error.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag {

template <typename coeff_t>
std::vector<int64_t> sparse_counts(kernels::CompiledOpSum const &ops,
                                   bool transpose) try {
  if (transpose) {
    // columns are shared between threads
    std::vector<int64_t> counts(size(ops.block_in), 0);
    kernels::elements<coeff_t>(ops, [&](int64_t, int64_t col, coeff_t) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++counts[col];
    });
    return counts;
  } else {
    // every row is counted by one thread
    std::vector<int64_t> counts(size(ops.block_out), 0);
    kernels::elements<coeff_t>(
        ops, [&](int64_t row, int64_t, coeff_t) { ++counts[row]; });
    return counts;
  }
}
XDIAG_CATCH

template <typename idx_t, typename coeff_t>
void sparse_fill(kernels::CompiledOpSum const &ops,
                 std::vector<int64_t> const &counts, idx_t i0, bool transpose,
                 idx_t *ptr, idx_t *idx, coeff_t *data) try {
  // offset[g]: the next free slot of group g, starting at its first entry
  int64_t ngroups = counts.size();
  std::vector<int64_t> offset(ngroups + 1, 0);
  std::partial_sum(counts.begin(), counts.end(), offset.begin() + 1);
  for (int64_t g = 0; g <= ngroups; ++g) {
    ptr[g] = (idx_t)(offset[g] + i0);
  }

  if (transpose) {
    // columns are shared between threads: claim the slots atomically
    kernels::elements<coeff_t>(ops, [&](int64_t row, int64_t col, coeff_t val) {
      int64_t k;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
      k = offset[col]++;
      idx[k] = (idx_t)(row + i0);
      data[k] = val;
    });
  } else {
    // every row is filled by one thread, in a fixed order
    kernels::elements<coeff_t>(ops, [&](int64_t row, int64_t col, coeff_t val) {
      int64_t k = offset[row]++;
      idx[k] = (idx_t)(col + i0);
      data[k] = val;
    });
  }

  // sort the entries of every group by idx (required for valid CSR / CSC)
  int64_t max_count =
      counts.empty() ? 0 : *std::max_element(counts.begin(), counts.end());
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    std::vector<int64_t> perm(max_count);
    std::vector<idx_t> idx_tmp(max_count);
    std::vector<coeff_t> data_tmp(max_count);
#ifdef _OPENMP
#pragma omp for
#endif
    for (int64_t g = 0; g < ngroups; ++g) {
      int64_t start = (int64_t)ptr[g] - i0;
      int64_t n = counts[g];
      std::copy(idx + start, idx + start + n, idx_tmp.begin());
      std::copy(data + start, data + start + n, data_tmp.begin());
      std::iota(perm.begin(), perm.begin() + n, (int64_t)0);
      std::sort(perm.begin(), perm.begin() + n,
                [&](int64_t a, int64_t b) { return idx_tmp[a] < idx_tmp[b]; });
      for (int64_t k = 0; k < n; ++k) {
        idx[start + k] = idx_tmp[perm[k]];
        data[start + k] = data_tmp[perm[k]];
      }
    }
  }
}
XDIAG_CATCH

template std::vector<int64_t>
sparse_counts<double>(kernels::CompiledOpSum const &, bool);
template std::vector<int64_t>
sparse_counts<complex>(kernels::CompiledOpSum const &, bool);

#define XDIAG_INST(IDX, COEFF)                                                 \
  template void sparse_fill<IDX, COEFF>(kernels::CompiledOpSum const &,        \
                                        std::vector<int64_t> const &, IDX,     \
                                        bool, IDX *, IDX *, COEFF *);
XDIAG_INST(int32_t, double)
XDIAG_INST(int32_t, complex)
XDIAG_INST(int64_t, double)
XDIAG_INST(int64_t, complex)
#undef XDIAG_INST

} // namespace xdiag
