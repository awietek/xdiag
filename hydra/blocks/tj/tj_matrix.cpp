#include "tj_matrix.h"

#include <hydra/blocks/tj/terms/apply_terms_dispatch.h>
#include <hydra/blocks/tj/terms/compile.h>
#include <hydra/operators/compiler.h>

namespace hydra {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, tJ const &block_in,
                              tJ const &block_out) {
  assert(block_in == block_out); // only temporary
  BondList bonds_c = tj::compile(bonds, 1e-12);
  operators::check_bonds_in_range(bonds, block_in.n_sites());

  if ((is_real<coeff_t>()) && (bonds_c.is_complex())) {
    Log.err("Error in matrix_gen: trying to create a real matrix from an "
            "intrisically complex BondList");
  }
  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  arma::Mat<coeff_t> mat(dim_out, dim_in, arma::fill::zeros);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  tj::apply_terms_dispatch<coeff_t>(bonds_c, indexing_in, indexing_out, fill);
  return mat;
}

template arma::Mat<double> matrix_gen<double>(BondList const &bonds,
                                              tJ const &block_in,
                                              tJ const &block_out);

template arma::Mat<complex> matrix_gen<complex>(BondList const &bonds,
                                                tJ const &block_in,
                                                tJ const &block_out);

arma::Mat<double> matrix_real(BondList const &bonds, tJ const &block_in,
                              tJ const &block_out) {
  return matrix_gen<double>(bonds, block_in, block_out);
}

arma::Mat<complex> matrix_cplx(BondList const &bonds, tJ const &block_in,
                               tJ const &block_out) {
  return matrix_gen<complex>(bonds, block_in, block_out);
}
  
} // namespace hydra
