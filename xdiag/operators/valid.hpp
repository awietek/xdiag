// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

// Validates an Op against the registry: known type, correct site count,
// disjoint sites (if required), matrix present (if required).
void check_valid(Op const &op);
void check_valid(Monomial const &mono);
void check_valid(OpSum const &ops);

// Lower-level helpers used in compilation and dispatch
void must_have_sites(Op const &op);
void must_not_have_sites(Op const &op);
void must_have_nsites(Op const &op, int64_t n);
void must_have_disjoint_sites(Op const &op);
void must_have_sites_in_range(Op const &op, int64_t l, int64_t u);
void must_have_matrix(Op const &op);
void must_not_have_matrix(Op const &op);

} // namespace xdiag::operators
