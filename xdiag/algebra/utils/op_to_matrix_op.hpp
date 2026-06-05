// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>

namespace xdiag::algebra {

// Converts a single Op to an equivalent "Matrix" Op carrying its explicit
// matrix. Used to turn named spin-1/2 operators (Sz, S+, S-, Sx, Sy, SzSz,
// Exchange, SdotS, ScalarChirality) into a uniform "Matrix" form for kernels
// that operate purely on matrices.
//
// Site/kron convention (mirrors combine_matrix_ops::embed_op):
//   For sites [s0, s1, ...], bit j of the local state index encodes the spin
//   at sites[j]. In armadillo kron(A, B): A acts on the high (outer) bit,
//   B on the low (inner) bit. So for sites [i, j]:  kron(j_op, i_op).
//
// TODO: this currently hardcodes spin-1/2 matrices. For bosons and generic
// spin-S the elementary matrices change with the local Hilbert space
// dimension `d`, so this function should eventually take `d` as a parameter
// and dispatch to the appropriate matrix representations.
Op op_to_matrix_op(Op const &op);

} // namespace xdiag::algebra
