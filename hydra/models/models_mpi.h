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

std::tuple<int, int, double>
get_exchange_s1_s2_Jhalf_mpi(Bond const &bond, Couplings const &couplings);

} // namespace hydra
