#pragma once

namespace hydra {

template <class bit_t> class Spinhalf;
template <class bit_t, class GroupAction> class SpinhalfSymmetric;
template <class bit_t> class SpinhalfMPI;

template <class bit_t> class tJ;
template <class bit_t, class GroupAction> class tJSymmetric;

template <class bit_t> class Electron;
template <class bit_t, class GroupAction> class ElectronSymmetric;

template <class Block> struct is_mpi_block_t {
  static constexpr bool value = false;
};

template <class Block>
constexpr bool is_mpi_block = is_mpi_block_t<Block>::value;
  
} // namespace hydra
