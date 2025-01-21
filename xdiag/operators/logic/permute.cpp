#include "permute.hpp"

#include <xdiag/operators/logic/isapprox.hpp>

namespace xdiag {

Op permute(Op const &op, Permutation const &perm) try {
  if (op.hassites()) {
    auto type = op.type();
    auto sites = op.sites();

    std::vector<int64_t> sites_permuted(sites.size());
    int64_t i = 0;
    for (int64_t s : sites) {
      if (s < 0) {
        XDIAG_THROW(fmt::format("Cannot permute Op: found a site which is "
                                "negative: \"{}\"\nOp:\n{}",
                                s, to_string(op)));
      } else if (s >= perm.size()) {
        XDIAG_THROW(
            fmt::format("Cannot permute Op: Op has a site with number \"{}\" "
                        "which exceeds the length of the Permutation \"{}\", "
                        "\nOp:\n{}\nPermutation:\n{}",
                        s, perm.size(), to_string(op), to_string(perm)));
      } else {
        sites_permuted[i++] = perm[s];
      }
    }

    if (op.hasmatrix()) {
      return Op(type, sites_permuted, op.matrix());
    } else {
      return Op(type, sites_permuted);
    }
  } else {
    return op;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum permute(OpSum const &ops, Permutation const &perm) try {
  OpSum ops_permuted;
  for (auto const &[cpl, op] : ops.plain()) {
    ops_permuted += cpl * permute(op, perm);
  }
  return ops_permuted;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
