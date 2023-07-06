//
// Created by Alex Wietek on 7/5/23.
//
#pragma once

#include <extern/armadillo/armadillo>
#include <hydra/common.h>

namespace hydra {

template <typename vec_t, typename mat_t>
void orthogonalize_inplace(vec_t &&v, mat_t const &Q, int64_t max_col = -1,
                           int iterations = 2);
// orthogonalize vector against columns of an orthonormal matrix
//
// Parameters:
// v:          vector to be orthogonalized, orthogonal vector on exit
// Q:          matrix to be orthogonalized against, assumed to be
//             orthonormal
// max_col:    columns 0,...max_col-1 are orthogonalized against,
//             (default: -1, all columns are orthogonalized agains)
// iterations: number of repeated orthogonalizations, default 2

template <typename vec_t, typename mat_t>
vec_t orthogonalize(vec_t const &v, mat_t const &Q, int64_t max_col = -1,
                    int iterations = 2);

// orthogonalize vector agains columns of an orthonormal matrix
//
// Parameters:
// v:          vector to be orthogonalized
// Q:          matrix to be orthogonalized against, assumed to be
//             orthonormal
// max_col:    columns 0,...max_col-1 are orthogonalized against,
//             (default: -1, all columns are orthogonalized agains)
// iterations: number of repeated orthogonalizations, default 2
//
// Returns:    orthogonal vector

} // namespace hydra
