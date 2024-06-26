#pragma once

#include <xdiag/operators/bond.hpp>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::operators {

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

// Standard checks of bond format
void check_bond_in_range(Bond const &bond, int64_t n_sites);
void check_bonds_in_range(BondList const &bonds, int64_t n_sites);
void check_bond_has_correct_number_of_sites(Bond const &bond, int64_t ns);
void check_bond_has_disjoint_sites(Bond const &bond);

} // namespace xdiag::operators
