#pragma once

#include <hydra/operators/bond.h>
#include <hydra/operators/bondlist.h>

namespace hydra::operators {

bool coupling_defined(Bond const &bond, BondList const &bonds);
bool matrix_defined(Bond const &bond, BondList const &bonds);

complex coupling(Bond const &bond, BondList const &bonds);
arma::cx_mat matrix(Bond const &bond, BondList const &bonds);

BondList compile_explicit_couplings(BondList const &bonds, double precision,
                                    std::string undefined_behavior);
BondList compile_explicit_matrices(BondList const &bonds, double precision,
                                   std::string undefined_behavior);
BondList compile_explicit(BondList const &bonds, double precision,
                          std::string undefined_behavior);

} // namespace hydra::operators
