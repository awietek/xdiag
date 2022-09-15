#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>

#include <hydra/blocks/spinhalf/spinhalf.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, arma::Col<coeff_t> const &vec_in,
           Spinhalf<bit_t> const &block_out, arma::Col<coeff_t> &vec_out);

} // namespace hydra
