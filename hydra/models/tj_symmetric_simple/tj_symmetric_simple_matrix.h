#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/tj_symmetric_simple/tj_symmetric_simple.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           tJSymmetricSimple<bit_t, GroupAction> const &block_in,
           tJSymmetricSimple<bit_t, GroupAction> const &block_out);

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           tJSymmetricSimple<bit_t, GroupAction> const &block_in,
           tJSymmetricSimple<bit_t, GroupAction> const &block_out);

} // namespace hydra
