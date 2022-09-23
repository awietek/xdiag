#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/tj/tj.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void apply(BondList const &bonds, tJ<bit_t> const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ<bit_t> const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
