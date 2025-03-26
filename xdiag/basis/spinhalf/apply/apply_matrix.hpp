#pragma once

#include <functional>

#include <xdiag/basis/spinhalf/apply/apply_term_diag.hpp>
#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_no_sym.hpp>
#include <xdiag/basis/spinhalf/apply/apply_term_offdiag_sym.hpp>

#include <xdiag/common.hpp>
#include <xdiag/operators/logic/non_branching_op.hpp>

namespace xdiag::basis::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class basis_t,
          class fill_f>
void apply_matrix(Coupling const &cpl, Op const &op, basis_t const &basis_in,
                  basis_t const &basis_out, fill_f fill) try {

  // Decompose into sum of non-branching operators
  auto ops_nb = operators::non_branching_ops<bit_t, coeff_t>(cpl, op);

  // Loop over sum of non-branching operators
  for (auto const &op_nb : ops_nb) {

    // Diagonal terms
    if (op_nb.isdiagonal()) {
      auto term_coeff = [&op_nb](bit_t spins) -> coeff_t {
        bit_t local_spins = op_nb.extract(spins);
        return op_nb.coeff(local_spins);
      };
      apply_term_diag<bit_t, coeff_t>(basis_in, term_coeff, fill);
    }

    // Offdiagonal terms
    else {
      auto non_zero_term = [&op_nb](bit_t spins) -> bool {
        bit_t local_spins = op_nb.extract(spins);
        return op_nb.non_zero_term(local_spins);
      };
      auto term_action = [&op_nb](bit_t spins) -> std::pair<bit_t, coeff_t> {
        bit_t local_spins = op_nb.extract(spins);
        auto [local_spins_new, coeff] = op_nb.state_coeff(local_spins);
        bit_t spins_new = op_nb.deposit(local_spins_new, spins);
        return {spins_new, coeff};
      };

      // Dispatch either symmetric of unsymmetric term application
      if constexpr (symmetric) {
        apply_term_offdiag_sym<bit_t, coeff_t>(
            basis_in, basis_out, non_zero_term, term_action, fill);
      } else {
        apply_term_offdiag_no_sym<bit_t, coeff_t>(
            basis_in, basis_out, non_zero_term, term_action, fill);
      }
    }

  } // for (auto const &op_nb : ops_nb)
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::spinhalf
