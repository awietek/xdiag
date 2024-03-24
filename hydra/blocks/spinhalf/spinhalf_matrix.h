#pragma once

#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/blocks/spinhalf/spinhalf.h>

namespace hydra {

arma::mat matrix(BondList const &bonds, Spinhalf const &block_in,
                 Spinhalf const &block_out);

arma::cx_mat matrixC(BondList const &bonds, Spinhalf const &block_in,
                     Spinhalf const &block_out);

// Developer functions
void matrix(double *mat, BondList const &bonds, Spinhalf const &block_in,
            Spinhalf const &block_out);

void matrixC(complex *mat, BondList const &bonds, Spinhalf const &block_in,
             Spinhalf const &block_out);

} // namespace hydra
