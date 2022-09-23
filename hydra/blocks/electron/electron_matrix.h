#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/electron/electron.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

template <typename bit_t, typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds,
                              Electron<bit_t> const &block_in,
                              Electron<bit_t> const &block_out);

template <typename bit_t>
inline arma::Mat<double> matrix_real(BondList const &bonds,
                                     Electron<bit_t> const &block_in,
                                     Electron<bit_t> const &block_out) {
  return matrix_gen<bit_t, double>(bonds, block_in, block_out);
}

template <typename bit_t>
inline arma::Mat<complex> matrix_cplx(BondList const &bonds,
                                      Electron<bit_t> const &block_in,
                                      Electron<bit_t> const &block_out) {
  return matrix_gen<bit_t, complex>(bonds, block_in, block_out);
}

} // namespace hydra
