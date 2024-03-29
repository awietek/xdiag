#pragma once

#include <hydra/extern/armadillo/armadillo>

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

#include <hydra/blocks/spinhalf/spinhalf.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, Spinhalf const &block_in,
           arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
