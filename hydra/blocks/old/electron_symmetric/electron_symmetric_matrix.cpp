#include "electron_symmetric_matrix.h"

#include <hydra/blocks/electron_symmetric/electron_symmetric.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_exchange.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_hopping.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_ising.h>
#include <hydra/blocks/electron_symmetric/terms/electron_symmetric_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <class bit_t, typename coeff_t>
lila::Matrix<coeff_t> MatrixGen(BondList const &bonds,
                                Couplings const &couplings,
                                ElectronSymmetric<bit_t> const &block_in,
                                ElectronSymmetric<bit_t> const &block_out) {

  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "electron_symmetric_matrix");
  idx_t dim = block_in.size();
  auto mat = lila::Zeros<coeff_t>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  using namespace terms;
  electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
  electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                             fill);
  electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &, Couplings const &,
                            ElectronSymmetric<uint16_t> const &,
                            ElectronSymmetric<uint16_t> const &);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &, Couplings const &,
                            ElectronSymmetric<uint32_t> const &,
                            ElectronSymmetric<uint32_t> const &);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &, Couplings const &,
                            ElectronSymmetric<uint64_t> const &,
                            ElectronSymmetric<uint64_t> const &);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &, Couplings const &,
                             ElectronSymmetric<uint16_t> const &,
                             ElectronSymmetric<uint16_t> const &);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &, Couplings const &,
                             ElectronSymmetric<uint32_t> const &,
                             ElectronSymmetric<uint32_t> const &);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &, Couplings const &,
                             ElectronSymmetric<uint64_t> const &,
                             ElectronSymmetric<uint64_t> const &);

} // namespace hydra
