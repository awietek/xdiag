#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/electron/electron.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, arma::Col<coeff_t> const &vec_in,
           Electron<bit_t> const &block_out, arma::Col<coeff_t> &vec_out);

} // namespace hydra
