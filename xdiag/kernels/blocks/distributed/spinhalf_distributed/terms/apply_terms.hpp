// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <xdiag/armadillo.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/mpi/buffer.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

#include <xdiag/kernels/blocks/distributed/spinhalf_distributed/terms/apply_exchange.hpp>
#include <xdiag/kernels/blocks/distributed/spinhalf_distributed/terms/apply_spsm.hpp>
#include <xdiag/kernels/blocks/distributed/spinhalf_distributed/terms/apply_sz.hpp>
#include <xdiag/kernels/blocks/distributed/spinhalf_distributed/terms/apply_szsz.hpp>
#include <xdiag/kernels/blocks/distributed/spinhalf_distributed/terms/transpose.hpp>

namespace xdiag::basis::spinhalf_distributed {

// Matrix-free, MPI-aware application of a (number/Sz conserving, no symmetry)
// OpSum on a distributed spin-1/2 basis. The ops are normal-ordered to single
// elementary operators (SzSz, Sz, Exchange, S+, S-, Id); two-site offdiagonal
// ops are classified by whether they act purely on the postfix bits (local),
// purely on the prefix bits (needs a transpose), or are mixed (needs a custom
// all-to-all). Anything beyond single-operator terms throws "not implemented".
template <typename coeff_t, class basis_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {

  auto algebra = algebra::spinhalf_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  int64_t n_postfix_bits = basis_in.n_postfix_bits();
  auto on_prefix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s >= n_postfix_bits; });
  };
  auto on_postfix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s < n_postfix_bits; });
  };

  std::vector<std::pair<Coeff, Op>> diagonal, prefix, postfix, mixed;
  for (auto const &[c, monomial] : ops_compiled) {
    if (monomial.size() != 1) {
      XDIAG_THROW("SpinhalfDistributed only supports single-operator terms "
                  "(monomials of length 1). Products of operators are not "
                  "implemented.");
    }
    Op op = monomial[0];
    std::string type = op.type();
    if ((type == "SzSz") || (type == "Sz") || (type == "Id")) {
      diagonal.push_back({c, op});
    } else if ((type == "Exchange") || (type == "S+") || (type == "S-")) {
      if (on_prefix(op)) {
        prefix.push_back({c, op});
      } else if (on_postfix(op)) {
        postfix.push_back({c, op});
      } else {
        mixed.push_back({c, op});
      }
    } else {
      XDIAG_THROW(fmt::format(
          "Unsupported Op type for SpinhalfDistributed block: \"{}\"", type));
    }
  }

  // Diagonal operators (purely local)
  for (auto const &[c, op] : diagonal) {
    std::string type = op.type();
    if (type == "SzSz") {
      apply_szsz(c, op, basis_in, vec_in, vec_out);
    } else if (type == "Sz") {
      apply_sz(c, op, basis_in, vec_in, vec_out);
    } else { // Id
      coeff_t cc = c.scalar().as<coeff_t>();
      for (int64_t idx = 0; idx < (int64_t)vec_in.size(); ++idx) {
        vec_out(idx) += cc * vec_in(idx);
      }
    }
  }

  // Postfix operators (local, off-diagonal)
  for (auto const &[c, op] : postfix) {
    std::string type = op.type();
    if (type == "Exchange") {
      apply_exchange_postfix(c, op, basis_in, vec_in, vec_out);
    } else { // S+ / S-
      apply_spsm_postfix(c, op, basis_in, vec_in, basis_out, vec_out);
    }
  }

  // Prefix operators: transpose to postfix|prefix order, act in the buffers,
  // transpose back, then accumulate into vec_out.
  if (!prefix.empty()) {
    int64_t buffer_size = std::max(basis_out.size_max(), basis_in.size_max());
    mpi::buffer.reserve<coeff_t>(buffer_size);
    coeff_t *send_buffer = mpi::buffer.send<coeff_t>();

    transpose(basis_in, vec_in.memptr(), false);
    for (auto const &[c, op] : prefix) {
      std::string type = op.type();
      if (type == "Exchange") {
        apply_exchange_prefix<basis_t, coeff_t>(c, op, basis_in);
      } else { // S+ / S-
        apply_spsm_prefix<basis_t, coeff_t>(c, op, basis_in, basis_out);
      }
    }
    transpose(basis_out, mpi::buffer.recv<coeff_t>(), true);
    for (int64_t idx = 0; idx < (int64_t)vec_out.size(); ++idx) {
      vec_out(idx) += send_buffer[idx];
    }
  }

  // Mixed operators (custom all-to-all)
  for (auto const &[c, op] : mixed) {
    std::string type = op.type();
    if (type == "Exchange") {
      apply_exchange_mixed(c, op, basis_in, vec_in, vec_out);
    } else {
      XDIAG_THROW(fmt::format("Unsupported mixed Op type for "
                              "SpinhalfDistributed block: \"{}\"",
                              type));
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::spinhalf_distributed
#endif
