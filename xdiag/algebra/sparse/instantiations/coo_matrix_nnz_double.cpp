// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/algebra/sparse/coo_matrix_nnz.cpp>

namespace xdiag::algebra {

#ifdef _OPENMP

template std::vector<int64_t> coo_matrix_nnz_thread<double>(OpSum const &,
                                                            Spinhalf const &,
                                                            Spinhalf const &);
template std::vector<int64_t>
coo_matrix_nnz_thread<double>(OpSum const &, tJ const &, tJ const &);
template std::vector<int64_t> coo_matrix_nnz_thread<double>(OpSum const &,
                                                            Electron const &,
                                                            Electron const &);
#else
template int64_t coo_matrix_nnz<double>(OpSum const &, Spinhalf const &,
                                        Spinhalf const &);
template int64_t coo_matrix_nnz<double>(OpSum const &, tJ const &, tJ const &);
template int64_t coo_matrix_nnz<double>(OpSum const &, Electron const &,
                                        Electron const &);

#endif
} // namespace xdiag::algebra
