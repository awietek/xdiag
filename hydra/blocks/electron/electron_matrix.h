#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/electron/electron.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::Mat<double> matrix_real(BondList const &bonds, Electron const &block_in,
                              Electron const &block_out);

arma::Mat<complex> matrix_cplx(BondList const &bonds, Electron const &block_in,
                               Electron const &block_out);

template <typename coeff_t>
arma::Mat<coeff_t> matrix_gen(BondList const &bonds, Electron const &block_in,
                              Electron const &block_out);

} // namespace hydra
