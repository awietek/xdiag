#pragma once

#include <hydra/common.h>
#include <hydra/blocks/blocks.h>
#include <tuple>

namespace hydra::detail {

template <class bit_t> class SpinhalfMPI;

template <> struct is_mpi_block_t<SpinhalfMPI<uint16>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint32>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint64>> {
  static constexpr bool value = true;
};

template <class bit_t> class ElectronMPI;

template <> struct is_mpi_block_t<ElectronMPI<uint16>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<ElectronMPI<uint32>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<ElectronMPI<uint64>> {
  static constexpr bool value = true;
};

} // namespace hydra::detail
