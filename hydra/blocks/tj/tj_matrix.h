#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/blocks/tj/tj.h>
#include <hydra/common.h>
#include <hydra/operators/bondlist.h>

namespace hydra {

arma::mat matrix(BondList const &bonds, tJ const &block_in,
                 tJ const &block_out);

arma::cx_mat matrixC(BondList const &bonds, tJ const &block_in,
                     tJ const &block_out);

// Developer functions
void matrix(double *mat, BondList const &bonds, tJ const &block_in,
            tJ const &block_out);

void matrixC(complex *mat, BondList const &bonds, tJ const &block_in,
             tJ const &block_out);

} // namespace hydra
