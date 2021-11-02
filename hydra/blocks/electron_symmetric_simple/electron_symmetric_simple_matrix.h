#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/blocks/electron_symmetric_simple/electron_symmetric_simple.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out);
  
template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_in,
           ElectronSymmetricSimple<bit_t, GroupAction> const &block_out);

} // namespace hydra
