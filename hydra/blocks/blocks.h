#pragma once

#include <variant>

#include <hydra/common.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/electron/electron.h>


// Forward declaration of blocks with default parameters
namespace hydra {

using Block = std::variant<Spinhalf, tJ, Electron>;

idx_t size(Block const &block);
idx_t hash(Block const &block);

// template <typename bit_t> class SpinhalfMPI;
// template <typename bit_t> class ElectronMPI;

} // namespace hydra

namespace hydra::mpi {
template <class Block> struct is_mpi_block_t {
  static constexpr bool value = false;
};

// template <> struct is_mpi_block_t<SpinhalfMPI<uint16_t>> {
//   static constexpr bool value = true;
// };
// template <> struct is_mpi_block_t<SpinhalfMPI<uint32_t>> {
//   static constexpr bool value = true;
// };
// template <> struct is_mpi_block_t<SpinhalfMPI<uint64_t>> {
//   static constexpr bool value = true;
// };

// template <> struct is_mpi_block_t<ElectronMPI<uint16_t>> {
//   static constexpr bool value = true;
// };
// template <> struct is_mpi_block_t<ElectronMPI<uint32_t>> {
//   static constexpr bool value = true;
// };
// template <> struct is_mpi_block_t<ElectronMPI<uint64_t>> {
//   static constexpr bool value = true;
// };

template <class Block>
constexpr bool is_mpi_block = is_mpi_block_t<Block>::value;

} // namespace hydra::mpi
