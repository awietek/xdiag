#pragma once

#include <xdiag/basis/tj_distributed/apply/apply_exchange.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_hopping.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_ising.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_number.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_raise_lower.hpp>
#include <xdiag/common.hpp>

namespace xdiag::basis::tj_distributed {

template <typename coeff_t, class BasisIn, class BasisOut>
void apply_terms(OpSum const &ops, BasisIn const &basis_in,
                 arma::Col<coeff_t> const &vec_in, BasisOut const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {
  (void)basis_out;

  using bit_t = typename BasisIn::bit_t;
  
  // Ops applied in up/dn order
  for (auto op : ops) {
    std::string type = op.type();
    if ((type == "ISING") || (type == "TJISING")) {
      tj_distributed::apply_ising<bit_t, coeff_t>(op, basis_in, vec_in.memptr(),
                                                  vec_out.memptr());
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      tj_distributed::apply_number<bit_t, coeff_t>(
          op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "EXCHANGE") {
      tj_distributed::apply_exchange<bit_t, coeff_t>(
          op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "HOPDN") {
      tj_distributed::apply_hopping<bit_t, coeff_t>(
          op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "HOPUP") {
      continue;
    } else {
      XDIAG_THROW(std::string("Unknown Op type for \"tJDistributed\" block: ") +
                  type);
    }
  }

  // Ops applied in dn/up order

  // Perform a transpose to dn/up order
  basis_in.transpose(vec_in.memptr());

  // after transpose transposed vector is stored in send_buffer of
  // mpi::buffer, hence we use this as new input vector
  coeff_t *vec_in_trans = mpi::buffer.send<coeff_t>();

  // the results of the application of terms is then written to the
  // mpi recv buffer
  coeff_t *vec_out_trans = mpi::buffer.recv<coeff_t>();

  for (auto op : ops) {
    std::string type = op.type();
    if (type == "HOPUP") {
      tj_distributed::apply_hopping<bit_t, coeff_t>(op, basis_in, vec_in_trans,
                                                    vec_out_trans);
    } else if ((type == "ISING") || (type == "TJISING") ||
               (type == "EXCHANGE") || (type == "HOPDN") ||
               (type == "NUMBERUP") || (type == "NUMBERDN")) {
      continue;
    } else {
      XDIAG_THROW(std::string("Unknown Op type for \"tJDistributed\" block: ") +
                  type);
    }
  }

  // Finally we transpose back to send_buffer ...
  basis_in.transpose_r(vec_out_trans);

  //  ... and fill the results to vec_out
  coeff_t *send = mpi::buffer.send<coeff_t>();
  for (int64_t i = 0; i < basis_out.size(); ++i) {
    vec_out[i] += send[i];
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

} // namespace xdiag::basis::tj_distributed
