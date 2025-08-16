#include "csr_matrix_fill.hpp"

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/utils/timing.hpp>

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
void csr_matrix_fill(OpSum const &ops, block_t const &block_in,
                     block_t const &block_out, int64_t nnz,
                     std::vector<idx_t> const &n_elements_in_row, idx_t *rowptr,
                     idx_t *col, coeff_t *data, idx_t i0, bool transpose) try {

  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);

  auto t0 = rightnow();
  int64_t max_entries_in_row =
      *std::max_element(n_elements_in_row.begin(), n_elements_in_row.end());

  arma::Col<idx_t> n_elements_in_row_offset(n_elements_in_row.size());
  std::exclusive_scan(n_elements_in_row.begin(), n_elements_in_row.end(),
                      n_elements_in_row_offset.begin(), 0);
  std::exclusive_scan(n_elements_in_row.begin(), n_elements_in_row.end(),
                      rowptr, 0);
  if (transpose) {
    rowptr[ncols] = nnz;
  } else {
    rowptr[nrows] = nnz;
  }
  timing(t0, rightnow(), "Sparse matrix generation: reductions", 1);

  t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  if (transpose) {
    auto fill = [&](int64_t c, int64_t r, coeff_t d) {
#ifdef _OPENMP
      idx_t idx;
#pragma omp atomic capture
      idx = n_elements_in_row_offset[c]++;
#else
      idx_t idx = n_elements_in_row_offset[c]++;
#endif
      col[idx] = (idx_t)r;
      data[idx] = d;
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  } else { // !transpose
    auto fill = [&](int64_t c, int64_t r, coeff_t d) {
#ifdef _OPENMP
      idx_t idx;
#pragma omp atomic capture
      idx = n_elements_in_row_offset[r]++;
#else
      idx_t idx = n_elements_in_row_offset[r]++;
#endif
      col[idx] = (idx_t)c;
      data[idx] = d;
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
  }
  timing(t0, rightnow(), "Sparse matrix generation: matrix assembly", 1);

  t0 = rightnow();
  int64_t nrowst = transpose ? ncols : nrows;

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
    std::vector<idx_t> indices(max_entries_in_row, 0);
    std::vector<idx_t> coltmp(max_entries_in_row, 0);
    std::vector<coeff_t> datatmp(max_entries_in_row, 0);

#ifdef _OPENMP
#pragma omp for
#endif
    for (idx_t row = 0; row < nrowst; ++row) {
      idx_t start = rowptr[row];
      idx_t end = rowptr[row + 1];
      idx_t nelems = end - start;

      std::copy(col + start, col + end, coltmp.begin());
      std::copy(data + start, data + end, datatmp.begin());
      std::iota(indices.begin(), indices.begin() + nelems, 0);
      std::sort(indices.begin(), indices.begin() + nelems,
                [&](idx_t i, idx_t j) { return coltmp[i] < coltmp[j]; });

      for (idx_t i = 0; i < nelems; ++i) {
        col[start + i] = coltmp[indices[i]];
        data[start + i] = datatmp[indices[i]];
      }
    }
#ifdef _OPENMP
  }
#endif
  timing(t0, rightnow(), "Sparse matrix generation: sorting", 1);

  if (i0 == 1) {
    auto t0 = rightnow();

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < nrowst + 1; ++i) {
      ++rowptr[i];
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < nnz; ++i) {
      ++col[i];
    }
    timing(t0, rightnow(), "Sparse matrix generation: adjusting 1-indexing", 1);
  }
}
XDIAG_CATCH

// int64_t, double
template void csr_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, std::vector<int64_t> const &, int64_t *,
                              int64_t *, double *, int64_t, bool);
template void csr_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              std::vector<int64_t> const &, int64_t *,
                              int64_t *, double *, int64_t, bool);
template void csr_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, std::vector<int64_t> const &, int64_t *,
                              int64_t *, double *, int64_t, bool);

// int64_t, complex
template void csr_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, std::vector<int64_t> const &, int64_t *,
                              int64_t *, complex *, int64_t, bool);
template void csr_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              std::vector<int64_t> const &, int64_t *,
                              int64_t *, complex *, int64_t, bool);
template void csr_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, std::vector<int64_t> const &, int64_t *,
                              int64_t *, complex *, int64_t, bool);

// int32_t, double
template void csr_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);

// int32_t, complex
template void csr_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, complex *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              std::vector<int32_t> const &, int32_t *,
                              int32_t *, complex *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, complex *, int32_t, bool);

} // namespace xdiag::algebra
