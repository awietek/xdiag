#include "spinhalf_matrix.h"

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/spinhalf/terms/apply_terms_dispatch.h>
#include <hydra/blocks/spinhalf/terms/compile_terms.h>

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

#include <hydra/operators/operator_qns.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> MatrixGen(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<bit_t> const &block_in,
                             Spinhalf<bit_t> const &block_out) {

  auto [bonds_c, couplings_c] =
      terms::spinhalf::compile_terms(bonds, couplings);

  int n_up_out = utils::spinhalf_nup(bonds_c, couplings_c, block_in);
  if (n_up_out != block_out.n_up())
    Log.err("Incompatible n_up in MatrixGen: {} != {}", n_up_out,
            block_out.n_up());

  utils::check_operator_works_with<coeff_t>(bonds_c, couplings_c,
                                            "spinhalf_matrix");
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  arma::Mat<coeff_t> mat(dim_out, dim_in, arma::fill::zeros);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  terms::spinhalf::apply_terms_dispatch<bit_t, coeff_t>(
      bonds_c, couplings_c, indexing_in, indexing_out, fill);
  return mat;
}

template arma::Mat<double>
MatrixGen<uint16_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint16_t> const &block_in,
                            Spinhalf<uint16_t> const &block_out);
template arma::Mat<double>
MatrixGen<uint32_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint32_t> const &block_in,
                            Spinhalf<uint32_t> const &block_out);
template arma::Mat<double>
MatrixGen<uint64_t, double>(BondList const &bonds, Couplings const &couplings,
                            Spinhalf<uint64_t> const &block_in,
                            Spinhalf<uint64_t> const &block_out);

template arma::Mat<complex>
MatrixGen<uint16_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint16_t> const &block_in,
                             Spinhalf<uint16_t> const &block_out);
template arma::Mat<complex>
MatrixGen<uint32_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint32_t> const &block_in,
                             Spinhalf<uint32_t> const &block_out);
template arma::Mat<complex>
MatrixGen<uint64_t, complex>(BondList const &bonds, Couplings const &couplings,
                             Spinhalf<uint64_t> const &block_in,
                             Spinhalf<uint64_t> const &block_out);

} // namespace hydra
