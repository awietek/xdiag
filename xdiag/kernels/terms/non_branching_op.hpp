// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <utility>
#include <vector>

#include <xdiag/armadillo.hpp>
#include <xdiag/bits/get_set.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>

namespace xdiag::kernels {

template <typename bit_t, typename coeff_t> class NonBranchingOp {
public:
  NonBranchingOp(std::vector<int64_t> const &sites,
                 arma::Mat<coeff_t> const &mat, int64_t d,
                 double precision = 1e-12);

  inline bool isdiagonal() const { return diagonal_; }
  inline std::pair<int64_t, coeff_t> state_coeff(int64_t s) const {
    return hops_[s];
  }

  // Defined inline (hot: called per non-zero term in the matrix-vector
  // product). extract reads the d-ary digits at sites_ into a linear index;
  // deposit writes a linear index back out into the state at sites_.
  inline int64_t extract(bit_t state) const {
    int64_t local = 0;
    for (int64_t i = (int64_t)sites_.size() - 1; i >= 0; --i) {
      local = local * d_ + bits::get(state, sites_[i]);
    }
    return local;
  }
  inline bit_t deposit(int64_t local, bit_t state) const {
    for (int64_t i = 0; i < (int64_t)sites_.size(); ++i) {
      bits::set(state, sites_[i], local % d_);
      local /= d_;
    }
    return state;
  }

private:
  std::vector<int64_t> sites_;
  int64_t d_;
  bool diagonal_;
  std::vector<std::pair<int64_t, coeff_t>> hops_;
};

template <typename bit_t, typename coeff_t>
std::vector<NonBranchingOp<bit_t, coeff_t>>
non_branching_ops(Coeff const &c, Op const &op, int64_t d,
                  double precision = 1e-12);

} // namespace xdiag::kernels
