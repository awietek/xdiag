// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "kernels.hpp"

#include <numeric>

#include <xdiag/armadillo.hpp>
#include <xdiag/basis/basis_onthefly.hpp>
#include <xdiag/basis/basis_symmetric.hpp>
#include <xdiag/bits/bitset.hpp>
#include <xdiag/combinatorics/combinations/combinations.hpp>
#include <xdiag/combinatorics/combinations/lin_table.hpp>
#include <xdiag/combinatorics/subsets/subsets.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/matrices/fill_functions.hpp>
#include <xdiag/matrices/spinhalf/matrix_policy.hpp>
#include <xdiag/utils/error.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace xdiag::matrices {

template <typename matrix_policy_t, typename basis_t, typename mat_t>
void apply(OpSum const &ops, basis_t const &basis_in, mat_t const &mat_in,
           basis_t const &basis_out, mat_t &mat_out) try {
  using coeff_t = typename mat_t::elem_type;
  mat_out.zeros();
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_apply(mat_in, mat_out, idx_in, idx_out, val);
      });
}
XDIAG_CATCH

template <typename matrix_policy_t, typename coeff_t, typename basis_t>
void matrix(OpSum const &ops, basis_t const &basis_in, basis_t const &basis_out,
            coeff_t *mat) try {
  int64_t m = basis_out.size();
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_matrix(mat, m, idx_in, idx_out, val);
      });
}
XDIAG_CATCH

#ifdef _OPENMP

template <typename matrix_policy_t, typename coeff_t, typename basis_t>
std::vector<int64_t> coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out) try {
  omp_set_schedule(omp_sched_static, 0);
  int nthreads = omp_get_max_threads();
  std::vector<int64_t> nnz_thread(nthreads, 0);
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out, [&](int64_t, int64_t, coeff_t, int num_thread) {
        fill_coo_count(nnz_thread, num_thread);
      });
  return nnz_thread;
}
XDIAG_CATCH

template <typename matrix_policy_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out,
                     std::vector<int64_t> const &nnz_thread, idx_t *rows,
                     idx_t *cols, coeff_t *data, idx_t i0) try {
  int nthreads = (int)nnz_thread.size();
  std::vector<int64_t> offset(nthreads, 0);
  std::partial_sum(nnz_thread.begin(), nnz_thread.end() - 1,
                   offset.begin() + 1);
  omp_set_schedule(omp_sched_static, 0);
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val, int num_thread) {
        fill_coo(rows, cols, data, offset, idx_in, idx_out, val, i0,
                 num_thread);
      });
}
XDIAG_CATCH

#else // !_OPENMP

template <typename matrix_policy_t, typename coeff_t, typename basis_t>
int64_t coo_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                       basis_t const &basis_out) try {
  int64_t nnz = 0;
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t, int64_t, coeff_t) { fill_coo_count(nnz); });
  return nnz;
}
XDIAG_CATCH

template <typename matrix_policy_t, typename coeff_t, typename basis_t,
          typename idx_t>
void coo_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, idx_t *rows, idx_t *cols,
                     coeff_t *data, idx_t i0) try {
  int64_t k = 0;
  matrix_policy_t::template call<coeff_t>(
      ops, basis_in, basis_out,
      [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
        fill_coo(rows, cols, data, k, idx_in, idx_out, val, i0);
      });
}
XDIAG_CATCH

#endif // _OPENMP

template <typename matrix_policy_t, typename coeff_t, typename basis_t>
std::vector<int64_t> csr_matrix_nnz(OpSum const &ops, basis_t const &basis_in,
                                    basis_t const &basis_out,
                                    bool transpose) try {
  if (transpose) {
    std::vector<int64_t> n_elements(basis_in.size(), 0);
    matrix_policy_t::template call<coeff_t>(
        ops, basis_in, basis_out, [&](int64_t idx_in, int64_t, coeff_t) {
          fill_csr_count(n_elements, idx_in);
        });
    return n_elements;
  } else {
    std::vector<int64_t> n_elements(basis_out.size(), 0);
    matrix_policy_t::template call<coeff_t>(
        ops, basis_in, basis_out, [&](int64_t, int64_t idx_out, coeff_t) {
          fill_csr_count(n_elements, idx_out);
        });
    return n_elements;
  }
}
XDIAG_CATCH

template <typename matrix_policy_t, typename coeff_t, typename basis_t,
          typename idx_t>
void csr_matrix_fill(OpSum const &ops, basis_t const &basis_in,
                     basis_t const &basis_out, std::vector<int64_t> &offset,
                     idx_t *col, coeff_t *data, idx_t i0, bool transpose) try {
  if (transpose) {
    // CSC: key by column (idx_in), store row (idx_out)
    matrix_policy_t::template call<coeff_t>(
        ops, basis_in, basis_out,
        [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
          fill_csr(offset, col, data, idx_out, idx_in, val, i0);
        });
  } else {
    // CSR: key by row (idx_out), store column (idx_in)
    matrix_policy_t::template call<coeff_t>(
        ops, basis_in, basis_out,
        [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
          fill_csr(offset, col, data, idx_in, idx_out, val, i0);
        });
  }
}
XDIAG_CATCH

} // namespace xdiag::matrices

// ---------------------------------------------------------------------------
// Explicit instantiations
// ---------------------------------------------------------------------------

using namespace arma;
using namespace xdiag;
using namespace xdiag::basis;
using namespace xdiag::bits;
using namespace xdiag::combinatorics;
using namespace xdiag::matrices;

#define INSTANTIATE_APPLY(MATRIX_POLICY, BASIS, MAT)                           \
  template void xdiag::matrices::apply<MATRIX_POLICY, BASIS, MAT>(             \
      OpSum const &, BASIS const &, MAT const &, BASIS const &, MAT &);

#define INSTANTIATE_MATRIX(MATRIX_POLICY, BASIS, COEFF)                        \
  template void xdiag::matrices::matrix<MATRIX_POLICY, COEFF, BASIS>(          \
      OpSum const &, BASIS const &, BASIS const &, COEFF *);

#ifdef _OPENMP
#define INSTANTIATE_COO_NNZ(MATRIX_POLICY, BASIS, COEFF)                       \
  template std::vector<int64_t>                                                \
  xdiag::matrices::coo_matrix_nnz<MATRIX_POLICY, COEFF, BASIS>(                \
      OpSum const &, BASIS const &, BASIS const &);
#define INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, IDX, COEFF)                 \
  template void                                                                \
  xdiag::matrices::coo_matrix_fill<MATRIX_POLICY, COEFF, BASIS, IDX>(          \
      OpSum const &, BASIS const &, BASIS const &,                             \
      std::vector<int64_t> const &, IDX *, IDX *, COEFF *, IDX);
#else
#define INSTANTIATE_COO_NNZ(MATRIX_POLICY, BASIS, COEFF)                       \
  template int64_t                                                             \
  xdiag::matrices::coo_matrix_nnz<MATRIX_POLICY, COEFF, BASIS>(                \
      OpSum const &, BASIS const &, BASIS const &);
#define INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, IDX, COEFF)                 \
  template void                                                                \
  xdiag::matrices::coo_matrix_fill<MATRIX_POLICY, COEFF, BASIS, IDX>(          \
      OpSum const &, BASIS const &, BASIS const &, IDX *, IDX *, COEFF *,      \
      IDX);
#endif

#define INSTANTIATE_CSR_NNZ(MATRIX_POLICY, BASIS, COEFF)                       \
  template std::vector<int64_t>                                                \
  xdiag::matrices::csr_matrix_nnz<MATRIX_POLICY, COEFF, BASIS>(                \
      OpSum const &, BASIS const &, BASIS const &, bool);
#define INSTANTIATE_CSR_FILL(MATRIX_POLICY, BASIS, IDX, COEFF)                 \
  template void                                                                \
  xdiag::matrices::csr_matrix_fill<MATRIX_POLICY, COEFF, BASIS, IDX>(          \
      OpSum const &, BASIS const &, BASIS const &, std::vector<int64_t> &,     \
      IDX *, COEFF *, IDX, bool);

#define INSTANTIATE_KERNELS(MATRIX_POLICY, BASIS)                              \
  INSTANTIATE_APPLY(MATRIX_POLICY, BASIS, vec)                                 \
  INSTANTIATE_APPLY(MATRIX_POLICY, BASIS, cx_vec)                              \
  INSTANTIATE_APPLY(MATRIX_POLICY, BASIS, mat)                                 \
  INSTANTIATE_APPLY(MATRIX_POLICY, BASIS, cx_mat)                              \
  INSTANTIATE_MATRIX(MATRIX_POLICY, BASIS, double)                             \
  INSTANTIATE_MATRIX(MATRIX_POLICY, BASIS, complex)                            \
  INSTANTIATE_COO_NNZ(MATRIX_POLICY, BASIS, double)                            \
  INSTANTIATE_COO_NNZ(MATRIX_POLICY, BASIS, complex)                           \
  INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, int32_t, double)                  \
  INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, int32_t, complex)                 \
  INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, int64_t, double)                  \
  INSTANTIATE_COO_FILL(MATRIX_POLICY, BASIS, int64_t, complex)                 \
  INSTANTIATE_CSR_NNZ(MATRIX_POLICY, BASIS, double)                            \
  INSTANTIATE_CSR_NNZ(MATRIX_POLICY, BASIS, complex)                           \
  INSTANTIATE_CSR_FILL(MATRIX_POLICY, BASIS, int32_t, double)                  \
  INSTANTIATE_CSR_FILL(MATRIX_POLICY, BASIS, int32_t, complex)                 \
  INSTANTIATE_CSR_FILL(MATRIX_POLICY, BASIS, int64_t, double)                  \
  INSTANTIATE_CSR_FILL(MATRIX_POLICY, BASIS, int64_t, complex)

//
//
// BasisOnTheFly Instantiations
//
//

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_subsets_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisOnTheFly<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_subsets_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisOnTheFly<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_lintable_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisOnTheFly<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_lintable_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisOnTheFly<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_bitset_dynamic)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_bitset_static_2)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_bitset_static_4)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_onthefly_combinations_bitset_static_8)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisOnTheFly<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP

//
//
// BasisSymmetric Instantiations
//
//

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_subsets_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSymmetric<Subsets<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_subsets_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSymmetric<Subsets<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_lintable_uint32_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSymmetric<LinTable<uint32_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_lintable_uint64_t)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSymmetric<LinTable<uint64_t>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_bitset_dynamic)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<BitsetDynamic>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_bitset_static_2)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<BitsetStatic2>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_bitset_static_4)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<BitsetStatic4>>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_symmetric_combinations_bitset_static_8)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy,
                    BasisSymmetric<Combinations<BitsetStatic8>>)
// END_INSTANTIATION_GROUP

//
//
// BasisSublattice Instantiations
//
//

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint32_t_1)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice32<1>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint32_t_2)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice32<2>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint32_t_3)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice32<3>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint32_t_4)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice32<4>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint32_t_5)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice32<5>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint64_t_1)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice64<1>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint64_t_2)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice64<2>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint64_t_3)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice64<3>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint64_t_4)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice64<4>)
// END_INSTANTIATION_GROUP

// BEGIN_INSTANTIATION_GROUP(spinhalf_sublattice_uint64_t_5)
INSTANTIATE_KERNELS(spinhalf::MatrixPolicy, BasisSublattice64<5>)
// END_INSTANTIATION_GROUP
