#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>
#include <xdiag/blocks/spinhalf/spinhalf.h>

namespace xdiag {

arma::mat matrix(BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out);

arma::cx_mat matrixC(BondList const &bonds, Spinhalf const &block_in,
                     Spinhalf const &block_out);

// Developer functions
void matrix(double *mat, BondList const &bonds, Spinhalf const &block_in,
            Spinhalf const &block_out);

void matrixC(complex *mat, BondList const &bonds, Spinhalf const &block_in,
             Spinhalf const &block_out);

} // namespace xdiag
