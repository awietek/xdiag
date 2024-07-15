#pragma once

#include <xdiag/operators/bond.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::operators {

BondList non_branching_bonds(Bond const &bond, double precision = 1e-12);
BondList non_branching_bonds(BondList const &bonds, double precision = 1e-12);
bool is_non_branching_bond(Bond const &bond, double precision = 1e-12);

template <typename bit_t, typename coeff_t> class NonBranchingBond {
public:
  explicit NonBranchingBond(Bond const &bond, double precision = 1e-12);
  bool is_diagonal() const;
  coeff_t coeff(bit_t local_state) const;
  std::pair<bit_t, coeff_t> state_coeff(bit_t local_state) const;
  bool non_zero_term(bit_t local_state) const;
  bit_t extract_local_state(bit_t state) const;
  bit_t deposit_local_state(bit_t local_state, bit_t state) const;
  int64_t number_difference() const;

private:
  std::vector<int64_t> sites_;
  bit_t dim_;
  bit_t mask_;
  arma::cx_mat matrix_;
  std::vector<bool> non_zero_term_;
  std::vector<bit_t> state_applied_;
  std::vector<coeff_t> coeff_;
};

} // namespace xdiag::operators
