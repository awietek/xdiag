// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::matrices {

template <typename bit_t, typename coeff_t> class NonBranchingOp {
public:
  NonBranchingOp(std::vector<int64_t> const &sites,
                 arma::Mat<coeff_t> const &mat, double precision = 1e-12);

  inline bool isdiagonal() const { return diagonal_; }
  inline std::pair<int64_t, coeff_t> state_coeff(int64_t s) const {
    return hops_[s];
  }

  int64_t extract(bit_t state) const;
  bit_t deposit(int64_t local, bit_t state) const;

private:
  std::vector<int64_t> sites_;
  bool diagonal_;
  std::vector<std::pair<int64_t, coeff_t>> hops_;
};

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coeff const &c, Op const &op, double precision = 1e-12);

} // namespace xdiag::matrices
