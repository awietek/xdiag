#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/electron/electron.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, Electron<bit_t> const &block_out);

template <class bit_t>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           Electron<bit_t> const &block_in, Electron<bit_t> const &block_out);

} // namespace hydra
