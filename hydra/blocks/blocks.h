#pragma once

namespace hydra {

template <class bit_t> class Spinhalf;
template <class bit_t, class GroupAction> class SpinhalfSymmetric;
template <class bit_t> class SpinhalfMPI;

template <class bit_t> class tJ;
template <class bit_t, class GroupAction> class tJSymmetric;
template <class bit_t, class GroupAction> class tJSymmetricSimple;

template <class bit_t> class Electron;
template <class bit_t, class GroupAction> class ElectronSymmetric;
template <class bit_t, class GroupAction> class ElectronSymmetricSimple;
template <class bit_t> class ElectronMPI;

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
