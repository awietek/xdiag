#pragma once

#include <hydra/common.h>
#include <hydra/symmetries/permutation_group_action.h>
#include <hydra/symmetries/permutation_group_lookup.h>

// Forward declaration of blocks with default parameters
namespace hydra {

template <typename bit_t = std_bit_t> class Spinhalf;
template <typename bit_t = std_bit_t> class SpinhalfSymmetric;
template <typename bit_t = std_bit_t> class SpinhalfMPI;

template <typename bit_t = std_bit_t> class tJ;
template <typename bit_t = std_bit_t> class tJSymmetric;

template <typename bit_t = std_bit_t> class Electron;
template <typename bit_t = std_bit_t,
          class GroupAction = PermutationGroupLookup<bit_t>>
class ElectronSymmetric;

template <typename bit_t = std_bit_t,
          class GroupAction = PermutationGroupAction>
class tJSymmetricSimple;
template <typename bit_t = std_bit_t,
          class GroupAction = PermutationGroupAction>
class ElectronSymmetricSimple;

template <typename bit_t = std_bit_t> class ElectronMPI;

} // namespace hydra

namespace hydra::detail {
template <class Block> struct is_mpi_block_t {
  static constexpr bool value = false;
};

template <> struct is_mpi_block_t<SpinhalfMPI<uint16_t>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint32_t>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<SpinhalfMPI<uint64_t>> {
  static constexpr bool value = true;
};

template <> struct is_mpi_block_t<ElectronMPI<uint16_t>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<ElectronMPI<uint32_t>> {
  static constexpr bool value = true;
};
template <> struct is_mpi_block_t<ElectronMPI<uint64_t>> {
  static constexpr bool value = true;
};

template <class Block>
constexpr bool is_mpi_block = is_mpi_block_t<Block>::value;

} // namespace hydra::detail
