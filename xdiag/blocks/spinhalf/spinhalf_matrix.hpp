#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/common.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag {

arma::mat matrix(BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out, double zero_precision = 1e-12);

arma::cx_mat matrixC(BondList const &bonds, Spinhalf const &block_in,
                     Spinhalf const &block_out, double zero_precision = 1e-12);

// Developer functions
void matrix(double *mat, BondList const &bonds, Spinhalf const &block_in,
            Spinhalf const &block_out, double zero_precision = 1e-12);

void matrixC(complex *mat, BondList const &bonds, Spinhalf const &block_in,
             Spinhalf const &block_out, double zero_precision = 1e-12);

} // namespace xdiag
