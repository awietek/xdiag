// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// Hermitian conjugate of a single Op. Returns an OpSum (rather than an Op)
// because hc may introduce a phase (e.g. hc(ExchangeAsym) = -ExchangeAsym).
// Throws if the Op type is unknown (not in the operator registry).
XDIAG_API OpSum hc(Op const &op);
XDIAG_API OpSum hc(Monomial const &mono);
XDIAG_API OpSum hc(OpSum const &ops);

// adjoint is the name of the same involution in the *-algebra language (it is
// the name used for it in Julia). It returns exactly hc(.).
XDIAG_API OpSum adjoint(Op const &op);
XDIAG_API OpSum adjoint(Monomial const &mono);
XDIAG_API OpSum adjoint(OpSum const &ops);

// Transpose and complex conjugate of the operator's matrix in the Fock basis.
// These are the two primitives; hermitian conjugation is hc = conj(transpose(.)).
//
//   transpose(A * B) = transpose(B) * transpose(A)
//   transpose(c * A) = c * transpose(A)
//   conj(A * B)      = conj(A) * conj(B)
//   conj(c * A)      = conj(c) * conj(A)
//
// Like hc, each returns an OpSum since a phase may appear (e.g. conj(Sy) = -Sy).
// Throws if the Op type is unknown (not in the operator registry).
XDIAG_API OpSum transpose(Op const &op);
XDIAG_API OpSum transpose(Monomial const &mono);
XDIAG_API OpSum transpose(OpSum const &ops);
XDIAG_API OpSum conj(Op const &op);
XDIAG_API OpSum conj(Monomial const &mono);
XDIAG_API OpSum conj(OpSum const &ops);

} // namespace xdiag
