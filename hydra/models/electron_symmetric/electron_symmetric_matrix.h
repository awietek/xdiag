#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron_symmetric/electron_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out);

} // namespace hydra
