// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "coo_matrix_nnz.hpp"

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/compilation.hpp>

namespace xdiag::algebra {

#ifdef _OPENMP

template <typename coeff_t, typename block_t>
std::vector<int64_t> coo_matrix_nnz_thread(OpSum const &ops,
                                           block_t const &block_in,
                                           block_t const &block_out) try {

  auto t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  int nthreads = 0;
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
  std::vector<int64_t> nnz_thread(nthreads, 0);

  omp_set_schedule(omp_sched_static, 0);
  auto fill_count_nnz = [&](int64_t idx_in, int64_t idx_out, coeff_t val,
                            int thread_num) { ++nnz_thread[thread_num]; };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_nnz);
  timing(t0, rightnow(), "Sparse COO: resource estimation", 1);
  return nnz_thread;
}
XDIAG_CATCH

template std::vector<int64_t> coo_matrix_nnz_thread<double>(OpSum const &,
                                                            Spinhalf const &,
                                                            Spinhalf const &);
template std::vector<int64_t> coo_matrix_nnz_thread<complex>(OpSum const &,
                                                             Spinhalf const &,
                                                             Spinhalf const &);

template std::vector<int64_t>
coo_matrix_nnz_thread<double>(OpSum const &, tJ const &, tJ const &);
template std::vector<int64_t>
coo_matrix_nnz_thread<complex>(OpSum const &, tJ const &, tJ const &);

template std::vector<int64_t> coo_matrix_nnz_thread<double>(OpSum const &,
                                                            Electron const &,
                                                            Electron const &);
template std::vector<int64_t> coo_matrix_nnz_thread<complex>(OpSum const &,
                                                             Electron const &,
                                                             Electron const &);

#else

template <typename coeff_t, typename block_t>
int64_t coo_matrix_nnz(OpSum const &ops, block_t const &block_in,
                       block_t const &block_out) try {
  auto t0 = rightnow();
  OpSum opsc = operators::compile<block_t>(ops);
  int64_t nnz = 0;
  auto fill_count_nnz = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    ++nnz;
  };
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill_count_nnz);
  timing(t0, rightnow(), "Sparse COO: resource estimation", 1);
  return nnz;
}
XDIAG_CATCH

template int64_t coo_matrix_nnz<double>(OpSum const &, Spinhalf const &,
                                        Spinhalf const &);
template int64_t coo_matrix_nnz<complex>(OpSum const &, Spinhalf const &,
                                         Spinhalf const &);

template int64_t coo_matrix_nnz<double>(OpSum const &, tJ const &, tJ const &);
template int64_t coo_matrix_nnz<complex>(OpSum const &, tJ const &, tJ const &);

template int64_t coo_matrix_nnz<double>(OpSum const &, Electron const &,
                                        Electron const &);
template int64_t coo_matrix_nnz<complex>(OpSum const &, Electron const &,
                                         Electron const &);

#endif

} // namespace xdiag::algebra
