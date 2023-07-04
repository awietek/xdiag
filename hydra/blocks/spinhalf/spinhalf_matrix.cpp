#include "spinhalf_matrix.h"

#include <hydra/blocks/spinhalf/spinhalf.h>
#include <hydra/blocks/spinhalf/terms/apply_terms_dispatch.h>
#include <hydra/blocks/spinhalf/terms/compile.h>
#include <hydra/blocks/spinhalf/terms/qns.h>
#include <hydra/blocks/utils/block_utils.h>
#include <hydra/operators/compiler.h>
#include <hydra/utils/logger.h>
#include <hydra/utils/timing.h>

namespace hydra {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Spinhalf const &block_in,
                              Spinhalf const &block_out) {

  BondList bonds_c = spinhalf::compile(bonds, 1e-12);
  operators::check_bonds_in_range(bonds, block_in.n_sites());

  if ((is_real<coeff_t>()) && (bonds_c.is_complex())) {
    Log.err("Error in matrix_gen: trying to create a real matrix from an "
            "intrisically complex BondList");
  }

  // if (block_in.n_up() != undefined_qn) {
  //   int n_up_out = spinhalf::nup(bonds_c) + block_in.n_up();
  //   if (n_up_out != block_out.n_up()) {
  //     Log.err("Incompatible n_up in matrix_gen: {} != {}", n_up_out,
  //             block_out.n_up());
  //   }
  // }

  idx_t dim_in = block_in.size();
  idx_t dim_out = block_out.size();

  arma::Mat<coeff_t> mat(dim_out, dim_in, arma::fill::zeros);
  auto fill = [&mat](idx_t idx_out, idx_t idx_in, coeff_t val) {
    mat(idx_out, idx_in) += val;
  };

  auto const &indexing_in = block_in.indexing();
  auto const &indexing_out = block_out.indexing();
  spinhalf::apply_terms_dispatch<coeff_t>(bonds_c, indexing_in, indexing_out,
                                          fill);
 return mat;
}

template arma::Mat<double> matrix_gen<double>(BondList const &bonds,
                                              Spinhalf const &block_in,
                                              Spinhalf const &block_out);

template arma::Mat<complex> matrix_gen<complex>(BondList const &bonds,
                                                Spinhalf const &block_in,
                                                Spinhalf const &block_out);

arma::Mat<double> matrix_real(BondList const &bonds, Spinhalf const &block_in,
                              Spinhalf const &block_out) {
  return matrix_gen<double>(bonds, block_in, block_out);
}

arma::Mat<complex> matrix_cplx(BondList const &bonds, Spinhalf const &block_in,
                               Spinhalf const &block_out) {
  return matrix_gen<complex>(bonds, block_in, block_out);
}

} // namespace hydra
