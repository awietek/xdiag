#include "electron_matrix.h"

#include <hydra/blocks/electron/terms/electron_exchange.h>
#include <hydra/blocks/electron/terms/electron_hopping.h>
#include <hydra/blocks/electron/terms/electron_ising.h>
#include <hydra/blocks/electron/terms/electron_u.h>

#include <hydra/blocks/electron/terms/electron_symmetric_exchange.h>
#include <hydra/blocks/electron/terms/electron_symmetric_hopping.h>
#include <hydra/blocks/electron/terms/electron_symmetric_ising.h>
#include <hydra/blocks/electron/terms/electron_symmetric_u.h>

#include <hydra/blocks/utils/block_utils.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> MatrixGen(BondList const &bonds, Couplings const &couplings,
                             Electron<bit_t> const &block_in,
                             Electron<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary

  utils::check_operator_works_with<coeff_t>(bonds, couplings,
                                            "electron_matrix");
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  arma::Mat<coeff_t> mat(dim_out, dim_in, arma::fill::zeros);

  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  if (block_in.symmetric()) {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_sym_np();
      terms::electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings,
                                                        indexing_in, fill);
      terms::electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings,
                                                      indexing_in, fill);
      terms::electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings,
                                                         indexing_in, fill);
    } else {
      auto const &indexing_in = block_in.indexing_sym_no_np();
      terms::electron_symmetric_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_symmetric_hopping<bit_t, coeff_t>(bonds, couplings,
                                                        indexing_in, fill);
      terms::electron_symmetric_ising<bit_t, coeff_t>(bonds, couplings,
                                                      indexing_in, fill);
      terms::electron_symmetric_exchange<bit_t, coeff_t>(bonds, couplings,
                                                         indexing_in, fill);
    }

  } else {

    if (block_in.charge_conserved() && block_in.sz_conserved()) {
      auto const &indexing_in = block_in.indexing_np();
      terms::electron_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
      terms::electron_ising<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                            fill);
      terms::electron_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                               fill);
    } else {
      auto const &indexing_in = block_in.indexing_no_np();
      terms::electron_U<bit_t, coeff_t>(couplings, indexing_in, fill);
      terms::electron_hopping<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                              fill);
      terms::electron_ising<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                            fill);
      terms::electron_exchange<bit_t, coeff_t>(bonds, couplings, indexing_in,
                                               fill);
    }
  }

  return mat;
}

template arma::Mat<double>
MatrixGen<uint16_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint16_t> const &block_in,
                            Electron<uint16_t> const &block_out);
template arma::Mat<double>
MatrixGen<uint32_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint32_t> const &block_in,
                            Electron<uint32_t> const &block_out);
template arma::Mat<double>
MatrixGen<uint64_t, double>(BondList const &bonds, Couplings const &couplings,
                            Electron<uint64_t> const &block_in,
                            Electron<uint64_t> const &block_out);

template arma::Mat<complex>
MatrixGen<uint16_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint16_t> const &block_in,
                             Electron<uint16_t> const &block_out);
template arma::Mat<complex>
MatrixGen<uint32_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint32_t> const &block_in,
                             Electron<uint32_t> const &block_out);
template arma::Mat<complex>
MatrixGen<uint64_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Electron<uint64_t> const &block_in,
                             Electron<uint64_t> const &block_out);

} // namespace hydra
