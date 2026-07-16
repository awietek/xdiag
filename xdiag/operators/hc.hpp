// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {

// Hermitian conjugate of a single Op. Returns an OpSum (rather than an Op)
// because hc may introduce a phase (e.g. hc(ExchangeAsym) = -ExchangeAsym).
// Throws if the Op type has no defined hermitian-conjugate partner.
XDIAG_API OpSum hc(Op const &op);
XDIAG_API OpSum hc(Monomial const &mono);
XDIAG_API OpSum hc(OpSum const &ops);

} // namespace xdiag
