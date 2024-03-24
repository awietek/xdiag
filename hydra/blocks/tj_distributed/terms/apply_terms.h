#pragma once

#include <hydra/common.h>

#include <hydra/blocks/tj_distributed/terms/apply_exchange.h>
#include <hydra/blocks/tj_distributed/terms/apply_hopping.h>
#include <hydra/blocks/tj_distributed/terms/apply_ising.h>
#include <hydra/blocks/tj_distributed/terms/apply_number.h>
#include <hydra/blocks/tj_distributed/terms/apply_raise_lower.h>

namespace hydra::tj_distributed {

template <typename bit_t, typename coeff_t, class BasisIn, class BasisOut>
void apply_terms(BondList const &bonds, BasisIn const &basis_in,
                 arma::Col<coeff_t> const &vec_in, BasisOut const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {
  (void)basis_out;

  // Bonds applied in up/dn order
  for (auto bond : bonds) {
    std::string type = bond.type();
    if ((type == "ISING") || (type == "TJISING")) {
      tj_distributed::apply_ising<bit_t, coeff_t>(
          bond, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if ((type == "NUMBERUP") || (type == "NUMBERDN")) {
      tj_distributed::apply_number<bit_t, coeff_t>(
          bond, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "EXCHANGE") {
      tj_distributed::apply_exchange<bit_t, coeff_t>(
          bond, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "HOPDN") {
      tj_distributed::apply_hopping<bit_t, coeff_t>(
          bond, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "HOPUP") {
      continue;
    } else {
      HydraThrow(
          std::runtime_error,
          std::string("Unknown bond type for \"tJDistributed\" block: ") +
              type);
    }
  }

  // Bonds applied in dn/up order

  // Perform a transpose to dn/up order
  basis_in.transpose(vec_in.memptr());

  // after transpose transposed vector is stored in send_buffer of
  // mpi::buffer, hence we use this as new input vector
  coeff_t *vec_in_trans = mpi::buffer.send<coeff_t>();

  // the results of the application of terms is then written to the
  // mpi recv buffer
  coeff_t *vec_out_trans = mpi::buffer.recv<coeff_t>();

  for (auto bond : bonds) {
    std::string type = bond.type();
    if (type == "HOPUP") {
      tj_distributed::apply_hopping<bit_t, coeff_t>(
          bond, basis_in, vec_in_trans, vec_out_trans);
    } else if ((type == "ISING") || (type == "TJISING") ||
               (type == "EXCHANGE") || (type == "HOPDN") ||
	       (type == "NUMBERUP") || (type == "NUMBERDN")) {
      continue;
    } else {
      HydraThrow(
          std::runtime_error,
          std::string("Unknown bond type for \"tJDistributed\" block: ") +
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

} catch (...) {
  HydraRethrow("Unable to apply terms for \"tJDistributed\" block");
}

} // namespace hydra::tj_distributed
