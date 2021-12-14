#pragma once

#include <lila/all.h>

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
lila::Matrix<coeff_t> MatrixGen(BondList const &bonds,
                                Couplings const &couplings,
                                tJSymmetric<bit_t> const &block_in,
                                tJSymmetric<bit_t> const &block_out);

template <typename bit_t>
lila::Matrix<double> MatrixReal(BondList const &bonds,
                                Couplings const &couplings,
                                tJSymmetric<bit_t> const &block_in,
                                tJSymmetric<bit_t> const &block_out) {
  return MatrixGen<bit_t, double>(bonds, couplings, block_in, block_out);
}

template <typename bit_t>
lila::Matrix<complex> MatrixCplx(BondList const &bonds,
                                 Couplings const &couplings,
                                 tJSymmetric<bit_t> const &block_in,
                                 tJSymmetric<bit_t> const &block_out) {
  return MatrixGen<bit_t, complex>(bonds, couplings, block_in, block_out);
}

} // namespace hydra
