#include "spinhalf_matrix.h"

#include <hydra/algebra/generic_operator.h>
#include <hydra/blocks/spinhalf/compile.h>
#include <hydra/operators/compiler.h>

namespace hydra {

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Spinhalf const &block_in,
                              Spinhalf const &block_out) try {
  return generic_matrix<coeff_t>(bonds, block_in, block_out, spinhalf::compile);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
  return arma::Mat<coeff_t>();
}

arma::mat matrix(BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out) try {
  return matrix_gen<double>(bonds, block_in, block_out);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
  return arma::mat();
}

arma::cx_mat matrixC(BondList const &bonds, Spinhalf const &block_in,
                     Spinhalf const &block_out) try {
  return matrix_gen<complex>(bonds, block_in, block_out);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
  return arma::cx_mat();
}

template <typename coeff_t>
void matrix_gen(coeff_t *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) try {
  generic_matrix<coeff_t>(mat, bonds, block_in, block_out, spinhalf::compile);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
}

void matrix(double *mat, BondList const &bonds, Spinhalf const &block_in,
            Spinhalf const &block_out) try {
  return matrix_gen<double>(mat, bonds, block_in, block_out);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
}

void matrixC(complex *mat, BondList const &bonds, Spinhalf const &block_in,
             Spinhalf const &block_out) try {
  return matrix_gen<complex>(mat, bonds, block_in, block_out);
} catch (...) {
  HydraRethrow("Cannot create matrix from Spinhalf block");
}

} // namespace hydra
