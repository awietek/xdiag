#pragma once

#include <hydra/models/models.h>
#include <tuple>

namespace hydra {

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
} // namespace hydra
