#include "symmetrize.hpp"

#include <xdiag/operators/logic/permute.hpp>

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group) try {
  auto ops = OpSum({op});
  return symmetrize(ops, group);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(Op const &op, Representation const &irrep) try {
  auto ops = OpSum({op});
  return symmetrize(ops, irrep);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(OpSum const &ops, PermutationGroup const &group) try {
  return symmetrize(ops, Representation(group));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T>
static OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                        arma::Col<T> const &characters) try {
  // TODO: needs to be adapted for MATRIX type (really? why?)
  OpSum ops_sym;
  int64_t N_group = group.size();
  for (auto [cpl, op] : ops.plain()) {

    std::string type = op.type();

    // Create all symmetrized ops
    for (int64_t i = 0; i < N_group; ++i) {
      Permutation perm = group[i];
      complex bloch = characters(i);
      Op op_perm = permute(op, perm);
      Scalar cpl_sym(bloch * cpl.scalar().as<complex>() / (complex)N_group);
      ops_sym += cpl_sym * op_perm;
    }
  }
  return ops_sym;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

OpSum symmetrize(OpSum const &ops, Representation const &irrep) try {
  if (isreal(irrep)) {
    auto characters = irrep.characters().as<arma::vec>();
    return symmetrize(ops, irrep.group(), characters);
  } else {
    auto characters = irrep.characters().as<arma::cx_vec>();
    return symmetrize(ops, irrep.group(), characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag
