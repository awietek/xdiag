#include "spinhalf_symmetric_matrix.h"

#include <hydra/utils/bitops.h>

#include <hydra/models/spinhalf_symmetric/terms/spinhalf_symmetric_exchange.h>
#include <hydra/models/spinhalf_symmetric/terms/spinhalf_symmetric_ising.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_in,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_out) {
  using namespace terms::spinhalf_symmetric;
  assert(block_in == block_out); // only temporary

  utils::check_symmetric_operator_real(
      bonds, couplings, block_in.irrep(), block_out.irrep(),
      "construct real SpinhalfSymmetric matrix");

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<double>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, block_in, fill);
  return mat;
}

template lila::Matrix<double> MatrixReal<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_out);
template lila::Matrix<double> MatrixReal<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_out);
template lila::Matrix<double> MatrixReal<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_out);

template lila::Matrix<double>
MatrixReal<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out);
template lila::Matrix<double>
MatrixReal<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out);
template lila::Matrix<double>
MatrixReal<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out);

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_in,
           SpinhalfSymmetric<bit_t, GroupAction> const &block_out) {
  using namespace terms::spinhalf_symmetric;
  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, block_in, fill);
  return mat;
}

template lila::Matrix<complex> MatrixCplx<uint16, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint16, PermutationGroupAction> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint32, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint32, PermutationGroupAction> const &block_out);
template lila::Matrix<complex> MatrixCplx<uint64, PermutationGroupAction>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_in,
    SpinhalfSymmetric<uint64, PermutationGroupAction> const &block_out);

template lila::Matrix<complex>
MatrixCplx<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    SpinhalfSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    SpinhalfSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    SpinhalfSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out);

} // namespace hydra
