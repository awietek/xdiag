#pragma once

#include <hydra/extern/armadillo/armadillo>

#include <hydra/blocks/tj/tj.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, tJ const &block_in,
           arma::Col<coeff_t> const &vec_in, tJ const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
