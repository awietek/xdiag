#pragma once

#include <xdiag/common.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

XDIAG_API OpSum symmetrize(Op const &op, PermutationGroup const &group);
XDIAG_API OpSum symmetrize(Op const &op, PermutationGroup const &group,
                           Representation const &irrep);
XDIAG_API OpSum symmetrize(OpSum const &ops, PermutationGroup const &group);
XDIAG_API OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                           Representation const &irrep);

} // namespace xdiag
