#pragma once

#include <functional>

#include <xdiag/blocks/spinhalf/terms/apply_term_offdiag_no_sym.hpp>
#include <xdiag/blocks/spinhalf/terms/apply_term_offdiag_sym.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/non_branching_bonds.hpp>

namespace xdiag::spinhalf {

template <typename bit_t, typename coeff_t, bool symmetric, class BasisIn,
          class BasisOut, class Fill>
void apply_non_branching(Bond const &bond, BasisIn &&basis_in,
                         BasisOut &&basis_out, Fill &&fill) {

  assert(bond.matrix_defined());
  assert(operators::is_non_branching_bond(bond));
  assert(bond.coupling_defined());

  auto bond_nb = operators::NonBranchingBond<bit_t, coeff_t>(bond);

  if (bond_nb.is_diagonal()) {
    auto term_coeff = [&bond_nb](bit_t spins) -> coeff_t {
      bit_t local_spins = bond_nb.extract_local_state(spins);
      return bond_nb.coeff(local_spins);
    };
    spinhalf::apply_term_diag<bit_t, coeff_t>(basis_in, term_coeff, fill);
  } else {
    auto non_zero_term = [&bond_nb](bit_t spins) -> bool {
      bit_t local_spins = bond_nb.extract_local_state(spins);
      return bond_nb.non_zero_term(local_spins);
    };
    auto term_action = [&bond_nb](bit_t spins) -> std::pair<bit_t, coeff_t> {
      bit_t local_spins = bond_nb.extract_local_state(spins);
      auto [local_spins_new, coeff] = bond_nb.state_coeff(local_spins);
      bit_t spins_new = bond_nb.deposit_local_state(local_spins_new, spins);

      // int64_t n_sites = 2;
      // Log("spins: {}, local_spins: {}, local_spins_new: {}, spins_new: {}",
      //     BSTR(spins), BSTR(local_spins), BSTR(local_spins_new),
      //     BSTR(spins_new));
      // std::cout << "coeff: "<< coeff<< std::endl;

      return {spins_new, coeff};
    };

    // Dispatch either symmetric of unsymmetric term application
    if constexpr (symmetric) {
      spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    } else {
      spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
          basis_in, basis_out, non_zero_term, term_action, fill);
    }
  }
}

} // namespace xdiag::spinhalf
