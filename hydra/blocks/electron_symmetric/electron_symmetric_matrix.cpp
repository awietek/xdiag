#include "electron_symmetric_matrix.h"

#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_hopping.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  using namespace terms::electron_symmetric;

  assert(block_in == block_out); // only temporary

  utils::check_symmetric_operator_real(
      bonds, couplings, block_in.irrep(), block_out.irrep(),
      "construct real ElectronSymmetric matrix");

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<double>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  do_U_symmetric(couplings, block_in, fill);
  do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16_t, PermutationGroupLookup<uint16_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16_t, PermutationGroupLookup<uint16_t>> const
        &block_in,
    ElectronSymmetric<uint16_t, PermutationGroupLookup<uint16_t>> const
        &block_out);
template lila::Matrix<double>
MatrixReal<uint32_t, PermutationGroupLookup<uint32_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32_t, PermutationGroupLookup<uint32_t>> const
        &block_in,
    ElectronSymmetric<uint32_t, PermutationGroupLookup<uint32_t>> const
        &block_out);
template lila::Matrix<double>
MatrixReal<uint64_t, PermutationGroupLookup<uint64_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64_t, PermutationGroupLookup<uint64_t>> const
        &block_in,
    ElectronSymmetric<uint64_t, PermutationGroupLookup<uint64_t>> const
        &block_out);

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  using namespace terms::electron_symmetric;

  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  do_U_symmetric(couplings, block_in, fill);
  do_hopping_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  return mat;
}

template lila::Matrix<complex>
MatrixCplx<uint16_t, PermutationGroupLookup<uint16_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16_t, PermutationGroupLookup<uint16_t>> const
        &block_in,
    ElectronSymmetric<uint16_t, PermutationGroupLookup<uint16_t>> const
        &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32_t, PermutationGroupLookup<uint32_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32_t, PermutationGroupLookup<uint32_t>> const
        &block_in,
    ElectronSymmetric<uint32_t, PermutationGroupLookup<uint32_t>> const
        &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64_t, PermutationGroupLookup<uint64_t>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64_t, PermutationGroupLookup<uint64_t>> const
        &block_in,
    ElectronSymmetric<uint64_t, PermutationGroupLookup<uint64_t>> const
        &block_out);

} // namespace hydra
