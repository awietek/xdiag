#pragma once

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t> class Spinhalf;
template <class bit_t> class SpinhalfMPI;

template <class bit_t> class tJ;

template <class bit_t> class Electron;
template <class bit_t, class SymmetryGroup> class ElectronSymmetric;

template <class Block> struct is_mpi_block_t {
  static constexpr bool value = false;
};

template <class Block>
constexpr bool is_mpi_block = is_mpi_block_t<Block>::value;

bool coupling_is_zero(Bond const &bond, Couplings const &couplings);

} // namespace hydra
