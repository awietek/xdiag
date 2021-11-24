#include "electron_matrix.h"

#include <hydra/blocks/electron/terms/electron_exchange.h>
#include <hydra/blocks/electron/terms/electron_hopping.h>
#include <hydra/blocks/electron/terms/electron_ising.h>
#include <hydra/blocks/electron/terms/electron_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
lila::Matrix<coeff_t>
MatrixGen(BondList const &bonds, Couplings const &couplings,
          Electron<bit_t> const &block_in, Electron<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings,
                                            "electron_matrix");
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  auto mat = lila::Zeros<coeff_t>(dim_out, dim_in);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  // auto const &indexing_out = block_out.indexing();

  terms::electron_U<bit_t, coeff_t>(couplings, indexing_in, fill);
  terms::electron_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  terms::electron_ising<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  terms::electron_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in, fill);
  return mat;
}

template lila::Matrix<double>
MatrixGen<uint16_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint16_t> const &block_in,
                            Electron<uint16_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint32_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint32_t> const &block_in,
                            Electron<uint32_t> const &block_out);
template lila::Matrix<double>
MatrixGen<uint64_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint64_t> const &block_in,
                            Electron<uint64_t> const &block_out);

template lila::Matrix<complex>
MatrixGen<uint16_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint16_t> const &block_in,
                             Electron<uint16_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint32_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint32_t> const &block_in,
                             Electron<uint32_t> const &block_out);
template lila::Matrix<complex>
MatrixGen<uint64_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint64_t> const &block_in,
                             Electron<uint64_t> const &block_out);

} // namespace hydra
