#include "tj_symmetric_matrix.h"

#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_exchange.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_hopping.h>
#include <hydra/blocks/tj_symmetric/terms/tj_symmetric_ising.h>
#include <hydra/blocks/tj_symmetric/tj_symmetric.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t>
lila::Matrix<double> MatrixReal(BondList const &bonds,
                                Couplings const &couplings,
                                tJSymmetric<bit_t> const &block_in,
                                tJSymmetric<bit_t> const &block_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary

  utils::check_symmetric_operator_real(bonds, couplings, block_in.irrep(),
                                       block_out.irrep(),
                                       "construct real tJSymmetric matrix");

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<double>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, double val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, double>(bonds, couplings, indexing_in, fill);

  return mat;
}

template <class bit_t>
lila::Matrix<complex> MatrixCplx(BondList const &bonds,
                                 Couplings const &couplings,
                                 tJSymmetric<bit_t> const &block_in,
                                 tJSymmetric<bit_t> const &block_out) {
  using namespace terms::tj_symmetric;

  assert(block_in == block_out); // only temporary

  idx_t dim = block_in.size();
  auto mat = lila::Zeros<complex>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, complex val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_hopping_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);
  do_ising_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, complex>(bonds, couplings, indexing_in, fill);

  return mat;
}

template lila::Matrix<double>
MatrixReal<uint16_t>(BondList const &, Couplings const &,
                     tJSymmetric<uint16_t> const &,
                     tJSymmetric<uint16_t> const &);
template lila::Matrix<double> MatrixReal<uint32>(BondList const &,
                                                 Couplings const &,
                                                 tJSymmetric<uint32> const &,
                                                 tJSymmetric<uint32> const &);
template lila::Matrix<double> MatrixReal<uint64>(BondList const &,
                                                 Couplings const &,
                                                 tJSymmetric<uint64> const &,
                                                 tJSymmetric<uint64> const &);

template lila::Matrix<complex>
MatrixCplx<uint16_t>(BondList const &, Couplings const &,
                     tJSymmetric<uint16_t> const &,
                     tJSymmetric<uint16_t> const &);
template lila::Matrix<complex>
MatrixCplx<uint32_t>(BondList const &, Couplings const &,
                     tJSymmetric<uint32_t> const &,
                     tJSymmetric<uint32_t> const &);
template lila::Matrix<complex>
MatrixCplx<uint64_t>(BondList const &, Couplings const &,
                     tJSymmetric<uint64_t> const &,
                     tJSymmetric<uint64_t> const &);

} // namespace hydra
