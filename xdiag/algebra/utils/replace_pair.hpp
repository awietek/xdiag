// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Given a monomial and a position k of an adjacent pair (k, k+1) that has been
// matched by a rule, builds the full replacement OpSum:
//   prefix * repl_pair * suffix
// where prefix = mono[0..k-1], suffix = mono[k+2..end].
OpSum replace_pair(Monomial const &mono, int64_t k, OpSum const &repl_pair);

} // namespace xdiag::algebra
