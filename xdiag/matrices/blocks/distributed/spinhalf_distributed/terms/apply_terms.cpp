// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_terms.hpp"

#include <xdiag/basis/apply_identity.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_exchange.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_spsm.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_sz.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_szsz.hpp>

#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>
#include <xdiag/basis/spinhalf_distributed/transpose.hpp>
#include <xdiag/utils/logger.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename coeff_t, class basis_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {
  // Determine diagonal and offdiagonal bonds
  OpSum ops_diagonal;
  OpSum ops_offdiagonal;
  for (auto [cpl, op] : ops.plain()) {
    if ((op.type() == "SzSz") || (op.type() == "Sz") || (op.type() == "Id")) {
      ops_diagonal += cpl * op;
    } else {
      ops_offdiagonal += cpl * op;
    }
  }

  // Among the offdiagonal bonds, determine which are acting on prefix, postfix
  // or are mixed
  OpSum ops_prefix, ops_postfix, ops_mixed;

  int64_t n_postfix_bits = basis_in.n_postfix_bits();
  auto isprefix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s >= n_postfix_bits; });
  };
  auto ispostfix = [&](Op const &op) {
    return std::all_of(op.sites().begin(), op.sites().end(),
                       [&](int64_t s) { return s < n_postfix_bits; });
  };

  for (auto [cpl, op] : ops_offdiagonal) {
    if (isprefix(op)) {
      ops_prefix += cpl * op;
    } else if (ispostfix(op)) {
      ops_postfix += cpl * op;
    } else {
      ops_mixed += cpl * op;
    }
  }

  // Diagonal operators
  double time_start = MPI_Wtime();
  for (auto [cpl, op] : ops_diagonal) {
    std::string type = op.type();
    if (type == "SzSz") {
      apply_szsz(cpl, op, basis_in, vec_in, vec_out);
    } else if (type == "Sz") {
      apply_sz(cpl, op, basis_in, vec_in, vec_out);
    } else if (type == "Id") {
      apply_identity<coeff_t>(cpl, basis_in,
                              [&](int64_t idxin, int64_t idxout, coeff_t c) {
                                vec_out[idxout] += c * vec_in[idxin];
                              });
    } else {
      XDIAG_THROW(fmt::format(
          "Unknown Op for SpinhalfDistributed block: \"{}\"", type));
    }
  }
  double time_end = MPI_Wtime();
  Log(3, "  diagonal ops: {:.6f} secs", time_end - time_start);

  // Apply postfix operators
  time_start = MPI_Wtime();
  for (auto [cpl, op] : ops_postfix) {
    std::string type = op.type();
    if (type == "Exchange") {
      apply_exchange_postfix(cpl, op, basis_in, vec_in, vec_out);
    } else if ((type == "S+") || (type == "S-")) {
      apply_spsm_postfix(cpl, op, basis_in, vec_in, basis_out, vec_out);
    } else {
      XDIAG_THROW(fmt::format(
          "Unknown Op for SpinhalfDistributed block: \"{}\"", type));
    }
  }
  time_end = MPI_Wtime();
  Log(3, "  postfix ops : {:.6f} secs", time_end - time_start);

  // Apply prefix operators (result is computed frmo send_buffer and stored in
  // recv_buffer)
  int64_t buffer_size = std::max(basis_out.size_max(), basis_in.size_max());
  mpi::buffer.reserve<coeff_t>(buffer_size);
  coeff_t *send_buffer = mpi::buffer.send<coeff_t>();
  coeff_t *recv_buffer = mpi::buffer.recv<coeff_t>();

  // Transpose to postfix | prefix order
  time_start = MPI_Wtime();
  transpose(basis_in, vec_in.memptr(), false);
  time_end = MPI_Wtime();
  Log(3, "  transpose   : {:.6f} secs", time_end - time_start);

  time_start = MPI_Wtime();
  for (auto [cpl, op] : ops_prefix) {
    std::string type = op.type();
    if (type == "Exchange") {
      apply_exchange_prefix<basis_t, coeff_t>(cpl, op, basis_in);
    } else if ((type == "S+") || (type == "S-")) {
      apply_spsm_prefix<basis_t, coeff_t>(cpl, op, basis_in, basis_out);
    } else {
      XDIAG_THROW(fmt::format(
          "Unknown Op for SpinhalfDistributed block: \"{}\"", type));
    }
  }
  time_end = MPI_Wtime();
  Log(3, "  prefix ops  : {:.6f} secs", time_end - time_start);

  // Transpose back to prefix | postfix order
  time_start = MPI_Wtime();
  transpose(basis_out, recv_buffer, true);
  time_end = MPI_Wtime();
  Log(3, "  transpose r : {:.6f} secs", time_end - time_start);

  // Fill contents of send buffer into vec_out
  time_start = MPI_Wtime();
  for (int64_t idx = 0; idx < vec_out.size(); ++idx) {
    vec_out(idx) += send_buffer[idx];
  }
  time_end = MPI_Wtime();
  Log(3, "  fill        : {:.6f} secs", time_end - time_start);

  /////////////////////////////
  // apply mixed operators
  time_start = MPI_Wtime();
  for (auto [cpl, op] : ops_mixed) {
    std::string type = op.type();
    if (type == "Exchange") {
      apply_exchange_mixed(cpl, op, basis_in, vec_in, vec_out);
    } else {
      XDIAG_THROW(fmt::format(
          "Unknown Op for SpinhalfDistributed block: \"{}\"", type));
    }
  }
  time_end = MPI_Wtime();
  Log(3, "  mixed       : {:.6f} secs", time_end - time_start);
}
XDIAG_CATCH

template void apply_terms(OpSum const &, BasisSz<uint32_t> const &,
                          arma::Col<double> const &, BasisSz<uint32_t> const &,
                          arma::Col<double> &);
template void apply_terms(OpSum const &, BasisSz<uint32_t> const &,
                          arma::Col<complex> const &, BasisSz<uint32_t> const &,
                          arma::Col<complex> &);
template void apply_terms(OpSum const &, BasisSz<uint64_t> const &,
                          arma::Col<double> const &, BasisSz<uint64_t> const &,
                          arma::Col<double> &);
template void apply_terms(OpSum const &, BasisSz<uint64_t> const &,
                          arma::Col<complex> const &, BasisSz<uint64_t> const &,
                          arma::Col<complex> &);

} // namespace xdiag::basis::spinhalf_distributed
