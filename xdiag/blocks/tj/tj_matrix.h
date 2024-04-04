#pragma once

#include <xdiag/extern/armadillo/armadillo>

#include <xdiag/common.h>
#include <xdiag/operators/bondlist.h>

#include <xdiag/blocks/tj/tj.h>

namespace xdiag {

arma::mat matrix(BondList const &bonds, tJ const &block_in,
                 tJ const &block_out);

arma::cx_mat matrixC(BondList const &bonds, tJ const &block_in,
                     tJ const &block_out);

// Developer functions
void matrix(double *mat, BondList const &bonds, tJ const &block_in,
            tJ const &block_out);

void matrixC(complex *mat, BondList const &bonds, tJ const &block_in,
             tJ const &block_out);

} // namespace xdiag
