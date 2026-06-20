// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include <xdiag/matrices/blocks/distributed/apply_distributed.hpp>

#include <xdiag/matrices/blocks/distributed/electron_distributed/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/distributed/electron_distributed/terms/apply_terms.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::matrices {

template <typename coeff_t>
void apply_distributed(OpSum const &ops, ElectronDistributed const &block_in,
                       arma::Col<coeff_t> const &vec_in,
                       ElectronDistributed const &block_out,
                       arma::Col<coeff_t> &vec_out) try {
  dispatch_basis(block_in, block_out,
                 [&](auto const &basis_in, auto const &basis_out) {
                   basis::electron_distributed::apply_terms<coeff_t>(
                       ops, basis_in, vec_in, basis_out, vec_out);
                 });
}
XDIAG_CATCH

template void apply_distributed(OpSum const &, ElectronDistributed const &,
                                arma::Col<double> const &, ElectronDistributed const &,
                                arma::Col<double> &);
template void apply_distributed(OpSum const &, ElectronDistributed const &,
                                arma::Col<complex> const &, ElectronDistributed const &,
                                arma::Col<complex> &);

} // namespace xdiag::matrices
