#pragma once

#include <optional>

#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

std::optional<Representation> representation(OpSum const &ops,
                                             PermutationGroup const &group);

} // namespace xdiag
