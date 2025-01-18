#include "qns.hpp"

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {
std::optional<Representation> representation(OpSum const &ops,
                                             PermutationGroup const &group) {
  std::vector<Scalar> characters;
  for (auto const &perm : group) {
    OpSum opsp = permute(ops, perm);
    std::optional<Scalar> factor = isapprox_multiple(opsp, ops);
    if (factor) {
      characters.push_back(*factor);
    } else {
      return std::nullopt;
    }
  }
  return Representation(group, characters);
}
} // namespace xdiag
