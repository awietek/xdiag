#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/tj_distributed/tj_distributed.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename coeff_t>
void apply(BondList const &bonds, tJDistributed const &block_in,
           arma::Col<coeff_t> const &vec_in, tJDistributed const &block_out,
           arma::Col<coeff_t> &vec_out);

} // namespace hydra
