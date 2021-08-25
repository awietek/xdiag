#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/tj_symmetric/tj_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t, class GroupAction>
void Apply(BondList const &bonds, Couplings const &couplings,
           tJSymmetric<bit_t, GroupAction> const &block_in,
           lila::Vector<complex> const &vec_in,
           tJSymmetric<bit_t, GroupAction> const &block_out,
           lila::Vector<complex> &vec_out);

} // namespace hydra
