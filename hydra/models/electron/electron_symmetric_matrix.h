#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron/electron_symmetric.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

// template <class bit_t>
// lila::Matrix<double> matrix_real(BondList const &bonds,
//                                  Couplings const &couplings,
//                                  ElectronSymmetric<bit_t> const &block_in,
//                                  ElectronSymmetric<bit_t> const &block_out);

template <class bit_t, class SymmetryGroup>
lila::Matrix<complex>
matrix_cplx(BondList const &bonds, Couplings const &couplings,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_in,
            ElectronSymmetric<bit_t, SymmetryGroup> const &block_out);

} // namespace hydra
