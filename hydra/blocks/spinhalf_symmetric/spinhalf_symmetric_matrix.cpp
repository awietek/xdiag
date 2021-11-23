#include "spinhalf_symmetric_matrix.h"

#include <hydra/bitops/bitops.h>

#include <hydra/blocks/spinhalf_symmetric/spinhalf_symmetric.h>
#include <hydra/blocks/spinhalf_symmetric/terms/spinhalf_symmetric_exchange.h>
#include <hydra/blocks/spinhalf_symmetric/terms/spinhalf_symmetric_ising.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
lila::Matrix<coeff_t> MatrixGen(BondList const &bonds,
                                Couplings const &couplings,
                                SpinhalfSymmetric<bit_t> const &block_in,
                                SpinhalfSymmetric<bit_t> const &block_out) {
  using namespace terms::spinhalf_symmetric;
  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings, block_in.irrep(),
                                            block_out.irrep(),
                                            "spinhalf_symmetric_apply");
  idx_t dim = block_in.size();
  auto mat = lila::Zeros<coeff_t>(dim, dim);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  do_ising_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  do_exchange_symmetric<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &, Couplings const &,
                            SpinhalfSymmetric<uint16_t> const &,
                            SpinhalfSymmetric<uint16_t> const &);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &, Couplings const &,
                            SpinhalfSymmetric<uint32_t> const &,
                            SpinhalfSymmetric<uint32_t> const &);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &, Couplings const &,
                            SpinhalfSymmetric<uint64_t> const &,
                            SpinhalfSymmetric<uint64_t> const &);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &, Couplings const &,
                             SpinhalfSymmetric<uint16_t> const &,
                             SpinhalfSymmetric<uint16_t> const &);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &, Couplings const &,
                             SpinhalfSymmetric<uint32_t> const &,
                             SpinhalfSymmetric<uint32_t> const &);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &, Couplings const &,
                             SpinhalfSymmetric<uint64_t> const &,
                             SpinhalfSymmetric<uint64_t> const &);

} // namespace hydra
