#pragma once

#include <xdiag/blocks/electron/electron.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag {

arma::mat matrix(BondList const &bonds, Electron const &block_in,
                 Electron const &block_out, double zero_precision=1e-12);

arma::cx_mat matrixC(BondList const &bonds, Electron const &block_in,
                     Electron const &block_out, double zero_precision=1e-12);

// Developer functions
void matrix(double *mat, BondList const &bonds, Electron const &block_in,
            Electron const &block_out, double zero_precision=1e-12);

void matrixC(complex *mat, BondList const &bonds, Electron const &block_in,
             Electron const &block_out, double zero_precision=1e-12);

} // namespace xdiag
