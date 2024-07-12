#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/basis/spinhalf_distributed/basis_sz.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/bondlist.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <class basis_t, typename coeff_t>
void apply_terms(BondList const &bonds, basis_t const &block_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &block_out,
                 arma::Col<coeff_t> &vec_out);

}

#endif
