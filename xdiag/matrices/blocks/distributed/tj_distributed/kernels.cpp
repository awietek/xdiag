// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/matrices/blocks/distributed/apply_distributed.hpp>

#include <xdiag/matrices/blocks/distributed/tj_distributed/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/distributed/tj_distributed/terms/apply_terms.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices {

template <typename coeff_t>
void apply_distributed(OpSum const &ops, tJDistributed const &block_in,
                       arma::Col<coeff_t> const &vec_in,
                       tJDistributed const &block_out,
                       arma::Col<coeff_t> &vec_out) try {
  dispatch_basis(block_in, block_out,
                 [&](auto const &basis_in, auto const &basis_out) {
                   basis::tj_distributed::apply_terms<coeff_t>(
                       ops, basis_in, vec_in, basis_out, vec_out);
                 });
}
XDIAG_CATCH

template void apply_distributed(OpSum const &, tJDistributed const &,
                                arma::Col<double> const &, tJDistributed const &,
                                arma::Col<double> &);
template void apply_distributed(OpSum const &, tJDistributed const &,
                                arma::Col<complex> const &, tJDistributed const &,
                                arma::Col<complex> &);

} // namespace xdiag::matrices
