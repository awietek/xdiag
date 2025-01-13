#include "symmetrize.hpp"

#include <xdiag/operators/logic/permute.hpp>

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group) try {
  auto ops = OpSum({op});
  return symmetrize(ops, group);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(Op const &op, PermutationGroup const &group,
                 Representation const &irrep) try {
  auto ops = OpSum({op});
  return symmetrize(ops, group, irrep);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(OpSum const &ops, PermutationGroup const &group) try {
  return symmetrize(ops, group, trivial_representation(group.size()));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// TODO: needs to be adapted for MATRIX type
OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                 Representation const &irrep) try {
  OpSum ops_sym;
  int64_t N_group = group.size();
  for (auto [cpl, op] : ops.plain()) {

    std::string type = op.type();

    // Create all symmetrized ops
    for (int64_t i = 0; i < N_group; ++i) {
      Permutation perm = group[i];
      complex bloch = irrep.character(i);
      Op op_perm = permute(op, perm);
      Scalar cpl_sym(bloch * cpl.scalar().as<complex>() / (complex)N_group);
      ops_sym += cpl_sym * op_perm;
    }
  }
  return ops_sym;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
