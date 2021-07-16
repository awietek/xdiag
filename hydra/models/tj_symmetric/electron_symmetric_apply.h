#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class SymmetryGroup>
void Apply(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
           lila::Vector<complex> const &vec_in,
           ElectronSymmetric<bit_t, SymmetryGroup> const &block_out,
           lila::Vector<complex> &vec_out);

} // namespace hydra
