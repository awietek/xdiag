#pragma once

#include <extern/armadillo/armadillo>
#include <hydra/blocks/blocks.h>
#include <hydra/operators/bondlist.h>
#include <hydra/states/state.h>

#include <hydra/blocks/electron/electron_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/tj/tj_apply.h>
#include <hydra/blocks/tj_distributed/tj_distributed_apply.h>

namespace hydra {

void apply(BondList const &bonds, State const &v, State &w);
  
template <typename coeff_t>
void apply(BondList const &bonds, block_variant_t const &block_in,
           arma::Col<coeff_t> const &vec_in, block_variant_t const &block_out,
           arma::Col<coeff_t> &vec_out);


} // namespace hydra
