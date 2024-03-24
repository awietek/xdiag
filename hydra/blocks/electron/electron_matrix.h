#pragma once

#include <hydra/extern/armadillo/armadillo>

#include <hydra/common.h>

#include <hydra/blocks/electron/electron.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::mat matrix(BondList const &bonds, Electron const &block_in,
                 Electron const &block_out);

arma::cx_mat matrixC(BondList const &bonds, Electron const &block_in,
                     Electron const &block_out);

// Developer functions
void matrix(double *mat, BondList const &bonds, Electron const &block_in,
            Electron const &block_out);

void matrixC(complex *mat, BondList const &bonds, Electron const &block_in,
             Electron const &block_out);

} // namespace hydra
