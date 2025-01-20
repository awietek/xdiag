#pragma once

#include <optional>

#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag {

std::optional<int64_t> nup(Op const &op);
std::optional<int64_t> nup(OpSum const &ops);
std::optional<int64_t> ndn(Op const &op);
std::optional<int64_t> ndn(OpSum const &ops);
std::optional<Representation> representation(OpSum const &ops,
                                             PermutationGroup const &group);

} // namespace xdiag
