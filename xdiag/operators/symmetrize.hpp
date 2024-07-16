#pragma once

#include <utility>

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group);
OpSum symmetrize(Op const &op, PermutationGroup const &group,
                 Representation const &irrep);
OpSum symmetrize(OpSum const &ops, PermutationGroup const &group);
OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                 Representation const &irrep);

} // namespace xdiag
