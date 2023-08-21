#include "spinhalf_matrix.h"

#include <hydra/algebra/generic_operator.h>

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
  return generic_matrix<coeff_t>(bonds, block_in, block_out, spinhalf::compile);
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

template <typename coeff_t>
void matrix_gen(coeff_t *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
  generic_matrix<coeff_t>(mat, bonds, block_in, block_out, spinhalf::compile);
}

template void matrix_gen<double>(double *mat, BondList const &bonds,
                                 Spinhalf const &block_in,
                                 Spinhalf const &block_out);

template void matrix_gen<complex>(complex *mat, BondList const &bonds,
                                  Spinhalf const &block_in,
                                  Spinhalf const &block_out);

void matrix_real(double *mat, BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out) {
  return matrix_gen<double>(mat, bonds, block_in, block_out);
}

void matrix_cplx(complex *mat, BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out) {
  return matrix_gen<complex>(mat, bonds, block_in, block_out);
}

} // namespace hydra
