#pragma once

#include <lila/all.h>

#include <hydra/blocks/electron/electron.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, lila::Vector<coeff_t> const &vec_in,
           Electron<bit_t> const &block_out, lila::Vector<coeff_t> &vec_out);

} // namespace hydra