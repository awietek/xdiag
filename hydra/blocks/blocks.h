#pragma once

#include <variant>

#include <hydra/common.h>
#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/tj/tj.h>
#include <hydra/blocks/electron/electron.h>

namespace hydra {

using block_variant_t = std::variant<Spinhalf, tJ, Electron>;
  
class Block {
public:
  Block() = default;
  Block(Spinhalf const& variant);
  Block(tJ const& variant);
  Block(Electron const& variant);
  Block(block_variant_t const& variant);
  idx_t size() const;
  uint64_t hash() const;

  bool operator==(Block const &rhs) const;
  bool operator!=(Block const &rhs) const;
  
  // Developer
  block_variant_t & variant();
  block_variant_t const& variant() const;
private:
  block_variant_t variant_;
};
  
  
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
