// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/algebra/sparse/coo_matrix_fill.cpp>

namespace xdiag::algebra {

#ifdef _OPENMP
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              std::vector<int64_t> const &, int64_t, int32_t *,
                              int32_t *, double *, int32_t);
#else
template void coo_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, int32_t *, int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              int32_t *, int32_t *, double *, int32_t);
template void coo_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, int32_t *, int32_t *, double *, int32_t);
#endif
} // namespace xdiag::algebra
