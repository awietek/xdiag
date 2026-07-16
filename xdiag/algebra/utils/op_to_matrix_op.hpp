// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>

namespace xdiag::algebra {

// Converts a single Op to an equivalent "Matrix" Op carrying its explicit
// matrix, for a local Hilbert space dimension `d`. Used to turn named operators
// into a uniform "Matrix" form for kernels that operate purely on matrices.
//
// The elementary matrices are built from `d`:
//   - Spin-S operators (Sz, S+, S-, Sx, Sy, SzSz, Exchange, ExchangeAsym,
//     SdotS, ScalarChirality) use the spin-S representation with S = (d-1)/2.
//     Basis index i = 0..d-1 corresponds to m = S - i (index 0 = highest m).
//     For d = 2 this reduces to the usual spin-1/2 matrices.
//   - Bosonic ladder/number operators (A, Adag, N) act on a Fock space
//     truncated to occupations 0..d-1, with basis index i = occupation i.
//
// Site/kron convention (mirrors combine_matrix_ops::embed_op):
//   For sites [s0, s1, ...], digit j of the local state index (base d) encodes
//   the local state at sites[j]. In armadillo kron(A, B): A acts on the high
//   (outer) digit, B on the low (inner) digit. So for sites [i, j]: kron(j_op,
//   i_op).
Op op_to_matrix_op(Op const &op, int64_t d);

} // namespace xdiag::algebra
