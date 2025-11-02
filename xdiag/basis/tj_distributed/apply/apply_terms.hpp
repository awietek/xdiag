// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/basis/tj_distributed/apply/apply_exchange.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_hopping.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_number.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_number_number.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_raise_lower.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_szsz.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t, class basis_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {
  using bit_t = typename basis_t::bit_t;

  // Adjust MPI buffer size if necessary
  int64_t buffer_size_in = std::max(basis_in.size(), basis_in.size_transpose());
  int64_t buffer_size_out =
      std::max(basis_out.size(), basis_out.size_transpose());
  mpi::buffer.reserve<coeff_t>(std::max(buffer_size_in, buffer_size_out));

  // Ops applied in up/dn order
  double time_start = MPI_Wtime();
  for (auto [cpl, op] : ops) {
    std::string type = op.type();
    if ((type == "SzSz") || (type == "tJSzSz")) {
      tj_distributed::apply_szsz<coeff_t>(cpl, op, basis_in, vec_in.memptr(),
                                          vec_out.memptr());
    } else if ((type == "Nup") || (type == "Ndn")) {
      tj_distributed::apply_number<coeff_t>(cpl, op, basis_in, vec_in.memptr(),
                                            vec_out.memptr());
    } else if (type == "NtotNtot") {
      tj_distributed::apply_ntot_ntot<coeff_t>(
          cpl, op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "Exchange") {
      tj_distributed::apply_exchange<coeff_t>(
          cpl, op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "Hopdn") {
      tj_distributed::apply_hopping<coeff_t>(cpl, op, basis_in, vec_in.memptr(),
                                             vec_out.memptr());
    } else if ((type == "Cdagdn") || (type == "Cdn")) {
      tj_distributed::apply_raise_lower<coeff_t>(
          cpl, op, basis_in, vec_in.memptr(), basis_out, vec_out.memptr());
    } else if ((type == "Hopup") || (type == "Cdagup") || (type == "Cup")) {
      continue;
    } else {
      XDIAG_THROW(std::string("Unknown Op type for \"tJDistributed\" block: ") +
                  type);
    }
  }
  double time_end = MPI_Wtime();
  Log(3, "  up/dn order : {:.6f} secs", time_end - time_start);

  // Ops applied in dn/up order

  // Perform a transpose to dn/up order
  time_start = MPI_Wtime();
  basis_in.transpose(vec_in.memptr());
  time_end = MPI_Wtime();
  Log(3, "  transpose   : {:.6f} secs", time_end - time_start);

  // after transpose transposed vector is stored in send_buffer of
  // mpi::buffer, hence we use this as new input vector
  coeff_t *vec_in_trans = mpi::buffer.send<coeff_t>();

  // the results of the application of terms is then written to the
  // mpi recv buffer
  coeff_t *vec_out_trans = mpi::buffer.recv<coeff_t>();

  time_start = MPI_Wtime();
  for (auto [cpl, op] : ops) {
    std::string type = op.type();
    if (type == "Hopup") {
      tj_distributed::apply_hopping<coeff_t>(cpl, op, basis_in, vec_in_trans,
                                             vec_out_trans);
    } else if ((type == "Cdagup") || (type == "Cup")) {
      tj_distributed::apply_raise_lower<coeff_t>(
          cpl, op, basis_in, vec_in_trans, basis_out, vec_out_trans);
    } else if ((type == "SzSz") || (type == "tJSzSz") || (type == "Exchange") ||
               (type == "Hopdn") || (type == "Nup") || (type == "Ndn") ||
               (type == "NtotNtot") || (type == "Cdagdn") || (type == "Cdn")) {
      continue;
    } else {
      XDIAG_THROW(std::string("Unknown Op type for \"tJDistributed\" block: ") +
                  type);
    }
  }
  time_end = MPI_Wtime();
  Log(3, "  dn/up order : {:.6f} secs", time_end - time_start);

  // Finally we transpose back to send_buffer ...
  time_start = MPI_Wtime();
  basis_out.transpose_r(vec_out_trans);
  time_end = MPI_Wtime();
  Log(3, "  transpose r : {:.6f} secs", time_end - time_start);

  //  ... and fill the results to vec_out
  coeff_t *send = mpi::buffer.send<coeff_t>();

  time_start = MPI_Wtime();
  for (int64_t i = 0; i < basis_out.size(); ++i) {
    vec_out[i] += send[i];
  }
  time_end = MPI_Wtime();
  Log(3, "  fill        : {:.6f} secs", time_end - time_start);
}
XDIAG_CATCH

} // namespace xdiag::basis::tj_distributed
