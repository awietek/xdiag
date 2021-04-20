#pragma once

#include <hydra/models/spinhalf.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

#include <lila/all.h>

namespace hydra {

template <class coeff_t, class bit_t = std_bit_t>
void apply(BondList bondlist, Couplings couplings,
           Spinhalf<bit_t> const &block_in, lila::Vector<coeff_t> const &vec_in,
           Spinhalf<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out);

} // namespace hydra
