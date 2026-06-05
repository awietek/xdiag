// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "symmetrize.hpp"

#include <xdiag/math/scalar.hpp>
#include <xdiag/algebra/permute.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

OpSum symmetrize(Op const &op, PermutationGroup const &group) try {
  auto ops = OpSum({op});
  return symmetrize(ops, group);
}
XDIAG_CATCH

OpSum symmetrize(Op const &op, Representation const &irrep) try {
  auto ops = OpSum({op});
  return symmetrize(ops, irrep);
}
XDIAG_CATCH

OpSum symmetrize(OpSum const &ops, PermutationGroup const &group) try {
  return symmetrize(ops, Representation(group));
}
XDIAG_CATCH

template <typename T>
static OpSum symmetrize(OpSum const &ops, PermutationGroup const &group,
                        arma::Col<T> const &characters) try {
  OpSum ops_sym;
  int64_t N_group = group.size();
  for (auto const &[coeff, mono] : ops.plain()) {
    for (int64_t i = 0; i < N_group; ++i) {
      Scalar cpls = coeff.scalar();
      T bloch = characters(i);
      Monomial mono_perm;
      for (auto const &op : mono) {
        mono_perm *= permute(op, group[i]);
      }
      Scalar cpl_sym = cpls * bloch / (T)N_group;
      ops_sym += cpl_sym * mono_perm;
    }
  }
  return ops_sym;
}
XDIAG_CATCH

OpSum symmetrize(OpSum const &ops, Representation const &irrep) try {
  if (!irrep.is_permutation()) {
    XDIAG_THROW(
        "symmetrize requires a Representation acting via a SitePermutation");
  }
  // Project with the conjugate characters conj(chi(g)) -- the standard
  // projection operator P = (1/|G|) sum_g conj(chi(g)) g -- so that the
  // symmetrized OpSum transforms under irrep (and not its conjugate). For real
  // characters conj is a no-op.
  PermutationGroup group = irrep.group();
  Vector characters = irrep.characters();
  if (isreal(irrep) && isreal(ops)) {
    return symmetrize(ops, group, characters.as<arma::vec>());
  } else {
    return symmetrize(ops, group,
                      arma::cx_vec(arma::conj(characters.as<arma::cx_vec>())));
  }
}
XDIAG_CATCH

} // namespace xdiag
