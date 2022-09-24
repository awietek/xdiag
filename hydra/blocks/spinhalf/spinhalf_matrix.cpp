#include "spinhalf_matrix.h"

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/spinhalf/terms/apply_terms_dispatch.h>
#include <hydra/blocks/spinhalf/terms/compile.h>
#include <hydra/blocks/spinhalf/terms/qns.h>

#include <hydra/blocks/utils/block_utils.h>
#include <hydra/utils/logger.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds,
                              Spinhalf<bit_t> const &block_in,
                              Spinhalf<bit_t> const &block_out) {

  BondList bonds_c = spinhalf::compile(bonds, 1e-12);

  if ((is_real<coeff_t>()) && (bonds_c.is_complex())) {
    Log.err("Error in matrix_gen: trying to create a real matrix from an "
            "intrisically complex BondList");
  }

  if (block_in.n_up() != undefined_qn) {
    int n_up_out = spinhalf::nup(bonds_c) + block_in.n_up();
    if (n_up_out != block_out.n_up()) {
      Log.err("Incompatible n_up in matrix_gen: {} != {}", n_up_out,
              block_out.n_up());
    }
  }
  
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  arma::Mat<coeff_t> mat(dim_out, dim_in, arma::fill::zeros);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  spinhalf::apply_terms_dispatch<bit_t, coeff_t>(bonds_c, indexing_in,
                                                 indexing_out, fill);
  return mat;
}

template arma::Mat<double>
matrix_gen<uint16_t, double>(BondList const &bonds,
                             Spinhalf<uint16_t> const &block_in,
                             Spinhalf<uint16_t> const &block_out);
template arma::Mat<double>
matrix_gen<uint32_t, double>(BondList const &bonds,
                             Spinhalf<uint32_t> const &block_in,
                             Spinhalf<uint32_t> const &block_out);
template arma::Mat<double>
matrix_gen<uint64_t, double>(BondList const &bonds,
                             Spinhalf<uint64_t> const &block_in,
                             Spinhalf<uint64_t> const &block_out);

template arma::Mat<complex>
matrix_gen<uint16_t, complex>(BondList const &bonds,
                              Spinhalf<uint16_t> const &block_in,
                              Spinhalf<uint16_t> const &block_out);
template arma::Mat<complex>
matrix_gen<uint32_t, complex>(BondList const &bonds,
                              Spinhalf<uint32_t> const &block_in,
                              Spinhalf<uint32_t> const &block_out);
template arma::Mat<complex>
matrix_gen<uint64_t, complex>(BondList const &bonds,
                              Spinhalf<uint64_t> const &block_in,
                              Spinhalf<uint64_t> const &block_out);

} // namespace hydra
