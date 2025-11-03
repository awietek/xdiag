#include "csr_matrix_nnz.hpp"

#include <numeric>

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/compilation.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::algebra {

template <typename idx_t, typename coeff_t, typename block_t>
std::vector<idx_t> csr_matrix_nnz(OpSum const &ops, block_t const &block_in,
                                  block_t const &block_out, bool transpose) {
  idx_t nrows = (idx_t)size(block_out);
  idx_t ncols = (idx_t)size(block_in);

  // Allocate memory for n_elements_in_row
  auto t0 = rightnow();
  std::vector<idx_t> n_elements_in_row;
  try {
    n_elements_in_row.resize(transpose ? ncols : nrows, 0);
  } catch (...) {
    XDIAG_THROW("Cannot allocate auxilliary for sparse matrix");
  }
  timing(t0, rightnow(),
         "Sparse matrix generation: auxiliary memory allocation", 1);

  // Count the number of non-zero elements in each row
  t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  if (transpose) {
    auto fill_count_row = [&](int64_t c, int64_t r, coeff_t val) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++n_elements_in_row[c];
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_row);

  } else { // not transposed
    auto fill_count_row = [&](int64_t c, int64_t r, coeff_t val) {
#ifdef _OPENMP
#pragma omp atomic update
#endif
      ++n_elements_in_row[r];
    };
    algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_row);
  }
  timing(t0, rightnow(), "Sparse matrix generation: resource estimation", 1);
  return n_elements_in_row;
}

} // namespace xdiag::algebra
