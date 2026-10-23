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

// Internal, not part of the public API: transpose and complex conjugate of the
// operator's matrix in the Fock basis. They are the two primitives of
// hc = conj(transpose(.)) and are used by the kernels (which run transpose(ops)).
// Unlike hc they are basis dependent: on a symmetry-adapted block with complex
// characters they relate a sector to its conjugate sector (momentum k to -k)
// rather than acting within the sector, which is why they are not exposed.
//
//   transpose(A * B) = transpose(B) * transpose(A)
//   transpose(c * A) = c * transpose(A)
//   conj(A * B)      = conj(A) * conj(B)
//   conj(c * A)      = conj(c) * conj(A)
//
// Like hc, each returns an OpSum since a phase may appear (e.g. conj(Sy) = -Sy).
// Throws if the Op type is unknown (not in the operator registry).
OpSum transpose(Op const &op);
OpSum transpose(Monomial const &mono);
OpSum transpose(OpSum const &ops);
OpSum conj(Op const &op);
OpSum conj(Monomial const &mono);
OpSum conj(OpSum const &ops);

} // namespace xdiag
