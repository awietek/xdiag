#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

#include <hydra/blocks/tj/tj.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
