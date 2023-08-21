#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/blocks.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::Mat<double> matrix_real(BondList const &bonds, Spinhalf const &block_in,
                              Spinhalf const &block_out);

arma::Mat<complex> matrix_cplx(BondList const &bonds, Spinhalf const &block_in,
                               Spinhalf const &block_out);

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Spinhalf const &block_in,
                              Spinhalf const &block_out);

void matrix_real(double *mat, BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out);

void matrix_cplx(complex *mat, BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out);

template <typename coeff_t>
void matrix_gen(coeff_t *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out);

} // namespace hydra
