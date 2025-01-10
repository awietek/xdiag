#pragma once

#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::operators {

template <typename bit_t, typename coeff_t> class NonBranchingOp {
public:
  explicit NonBranchingOp(std::vector<int64_t> const &sites,
                          arma::Mat<coeff_t> const &mat,
                          double precision = 1e-12);
  bool isdiagonal() const;
  coeff_t coeff(bit_t local_state) const;
  std::pair<bit_t, coeff_t> state_coeff(bit_t local_state) const;
  bool non_zero_term(bit_t local_state) const;

  bit_t extract(bit_t state) const;
  bit_t deposit(bit_t local_state, bit_t state) const;

private:
  bit_t mask_;
  bool diagonal_;
  std::vector<bool> non_zero_term_;
  std::vector<bit_t> state_applied_;
  std::vector<coeff_t> coeff_;
};

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coupling const &cpl, Op const &op, double precision = 1e-12);

} // namespace xdiag::operators
