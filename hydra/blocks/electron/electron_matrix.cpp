#include "electron_matrix.h"

#include <hydra/blocks/electron/terms/apply_terms_dispatch.h>
#include <hydra/blocks/electron/terms/compile.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds,
                              Electron<bit_t> const &block_in,
                              Electron<bit_t> const &block_out) {
  assert(block_in == block_out); // only temporary

  BondList bonds_c = electron::compile(bonds, 1e-12);

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
  electron::apply_terms_dispatch<bit_t, coeff_t>(bonds_c, indexing_in,
                                                 indexing_out, fill);
  return mat;
}

template arma::Mat<double>
matrix_gen<uint16_t, double>(BondList const &bonds,
                             Electron<uint16_t> const &block_in,
                             Electron<uint16_t> const &block_out);
template arma::Mat<double>
matrix_gen<uint32_t, double>(BondList const &bonds,
                             Electron<uint32_t> const &block_in,
                             Electron<uint32_t> const &block_out);
template arma::Mat<double>
matrix_gen<uint64_t, double>(BondList const &bonds,
                             Electron<uint64_t> const &block_in,
                             Electron<uint64_t> const &block_out);

template arma::Mat<complex>
matrix_gen<uint16_t, complex>(BondList const &bonds,
                              Electron<uint16_t> const &block_in,
                              Electron<uint16_t> const &block_out);
template arma::Mat<complex>
matrix_gen<uint32_t, complex>(BondList const &bonds,
                              Electron<uint32_t> const &block_in,
                              Electron<uint32_t> const &block_out);
template arma::Mat<complex>
matrix_gen<uint64_t, complex>(BondList const &bonds,
                              Electron<uint64_t> const &block_in,
                              Electron<uint64_t> const &block_out);

} // namespace hydra
