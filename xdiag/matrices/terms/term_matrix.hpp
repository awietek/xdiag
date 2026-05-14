// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/matrices/terms/non_branching_op.hpp>
#include <xdiag/matrices/terms/term_diag.hpp>
#include <xdiag/matrices/terms/term_offdiag.hpp>

namespace xdiag::matrices {

template <typename coeff_t, class basis_t, class fill_f>
void term_matrix(Coeff const &c, Op const &op, basis_t const &basis_in,
                 basis_t const &basis_out, fill_f fill) {
  using bit_t = typename basis_t::bit_t;

  // Decompose into sum of non-branching operators
  auto ops_nb = non_branching_ops<bit_t, coeff_t>(c, op);

  // Loop over sum of non-branching operators
  for (NonBranchingOp<bit_t, coeff_t> const &op_nb : ops_nb) {

    // Diagonal terms
    if (op_nb.isdiagonal()) {
      auto term_coeff = [&](bit_t spins) -> coeff_t {
        int64_t local = op_nb.extract(spins);
        return op_nb.state_coeff(local).second;
      };
      term_diag(basis_in, basis_out, term_coeff, fill);
    } else { // Offdiagonal terms
      auto non_zero_term = [&](bit_t spins) -> bool {
        int64_t local = op_nb.extract(spins);
        return op_nb.state_coeff(local).second != coeff_t(0);
      };
      auto term_action = [&](bit_t spins) -> std::pair<bit_t, coeff_t> {
        int64_t local = op_nb.extract(spins);
        auto [local_new, coeff] = op_nb.state_coeff(local);
        bit_t spins_new = op_nb.deposit(local_new, spins);
        return {spins_new, coeff};
      };
      term_offdiag(basis_in, basis_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace xdiag::matrices
