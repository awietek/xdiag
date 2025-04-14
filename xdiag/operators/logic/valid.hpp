// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {
void check_valid(Op const &op);
void check_valid(OpSum const &ops);

void check_valid(Op const &op, int64_t nsites);
void check_valid(OpSum const &ops, int64_t nsites);

void must_have_sites(Op const &op);
void must_not_have_sites(Op const &op);

void must_have_nsites(Op const &op, int64_t n);
void must_have_disjoint_sites(Op const &op);
void must_have_sites_in_range(Op const &op, int64_t l, int64_t u);

void must_have_matrix(Op const &op);
void must_not_have_matrix(Op const &op);

} // namespace xdiag
