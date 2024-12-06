#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag {
void check_valid(Op const &op);
void check_valid(OpSum const &ops);

void must_have_sites(Op const &op);
void must_not_have_sites(Op const &op);

void must_have_n_sites(Op const &op, int64_t n);
void must_have_disjoint_sites(Op const &op);

void must_have_matrix(Op const &op);
void must_not_have_matrix(Op const &op);

} // namespace xdiag
