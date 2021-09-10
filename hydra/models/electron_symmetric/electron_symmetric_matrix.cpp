#include "electron_symmetric_matrix.h"

// #include
// <hydra/models/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_hopping.h>
// #include <hydra/models/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/models/electron_symmetric/terms/electron_symmetric_u.h>

#include <hydra/models/utils/model_utils.h>

namespace hydra {

template <class bit_t, class GroupAction>
lila::Matrix<double>
MatrixReal(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  assert(block_in == block_out); // only temporary

  utils::check_symmetric_operator_real(
      bonds, couplings, block_in.irrep(), block_out.irrep(),
      "construct real ElectronSymmetric matrix");

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<double>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  electronterms::do_U_symmetric(couplings, block_in, fill);
  electronterms::do_hopping_symmetric<bit_t, double>(bonds, couplings, block_in,
                                                     fill);
  // electronterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  // electronterms::do_exchange_symmetric<bit_t, double>(bonds, couplings,
  // block_in,
  //                                                fill);
  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out);
template lila::Matrix<double>
MatrixReal<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out);
template lila::Matrix<double>
MatrixReal<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out);

template <class bit_t, class GroupAction>
lila::Matrix<complex>
MatrixCplx(BondList const &bonds, Couplings const &couplings,
           ElectronSymmetric<bit_t, GroupAction> const &block_in,
           ElectronSymmetric<bit_t, GroupAction> const &block_out) {
  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  electronterms::do_U_symmetric(couplings, block_in, fill);
  electronterms::do_hopping_symmetric<bit_t, complex>(bonds, couplings,
                                                      block_in, fill);
  // electronterms::do_ising_symmetric<bit_t>(bonds, couplings, block_in, fill);
  // electronterms::do_exchange_symmetric<bit_t, complex>(bonds, couplings,
  // block_in,
  //                                                 fill);
  return mat;
}

template lila::Matrix<complex>
MatrixCplx<uint16, PermutationGroupLookup<uint16>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_in,
    ElectronSymmetric<uint16, PermutationGroupLookup<uint16>> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint32, PermutationGroupLookup<uint32>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_in,
    ElectronSymmetric<uint32, PermutationGroupLookup<uint32>> const &block_out);
template lila::Matrix<complex>
MatrixCplx<uint64, PermutationGroupLookup<uint64>>(
    BondList const &bonds, Couplings const &couplings,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_in,
    ElectronSymmetric<uint64, PermutationGroupLookup<uint64>> const &block_out);

} // namespace hydra
