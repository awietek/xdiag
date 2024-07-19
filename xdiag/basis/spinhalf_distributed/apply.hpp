#pragma once
#ifdef XDIAG_USE_MPI

#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::basis::spinhalf_distributed {

template <typename coeff_t, class basis_t>
void apply(OpSum const &ops, basis_t const &basis_in,
           arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
           arma::Col<coeff_t> &vec_out);

}
#endif
