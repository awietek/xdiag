#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/blocks/electron_symmetric_simple/electron_symmetric_simple.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<double> const &vec_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<double> &vec_out);

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out);

} // namespace hydra
