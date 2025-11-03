// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/algebra/sparse/csr_matrix_fill.cpp>

namespace xdiag::algebra {
template void csr_matrix_fill(OpSum const &, Spinhalf const &, Spinhalf const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, tJ const &, tJ const &, int64_t,
                              std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);
template void csr_matrix_fill(OpSum const &, Electron const &, Electron const &,
                              int64_t, std::vector<int32_t> const &, int32_t *,
                              int32_t *, double *, int32_t, bool);
} // namespace xdiag::algebra
