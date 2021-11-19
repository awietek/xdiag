#pragma once

#include <lila/all.h>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double> MatrixReal(BondList const &bonds,
                                Couplings const &couplings,
                                tJSymmetric<bit_t> const &block_in,
                                tJSymmetric<bit_t> const &block_out);

template <class bit_t>
lila::Matrix<complex> MatrixCplx(BondList const &bonds,
                                 Couplings const &couplings,
                                 tJSymmetric<bit_t> const &block_in,
                                 tJSymmetric<bit_t> const &block_out);

} // namespace hydra
