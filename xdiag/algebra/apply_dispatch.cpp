// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "apply_dispatch.hpp"

#ifdef XDIAG_USE_MPI
#include <xdiag/basis/electron_distributed/apply/apply_terms.hpp>
#include <xdiag/basis/spinhalf_distributed/apply/apply_terms.hpp>
#include <xdiag/basis/tj_distributed/apply/apply_terms.hpp>
#endif

namespace xdiag::algebra {

#ifdef XDIAG_USE_MPI
template <typename coeff_t>
void apply_dispatch(OpSum const &ops, ElectronDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    ElectronDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  using namespace basis::electron_distributed;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      [&](auto &&basis_in, auto &&basis_out) {
        using basis_in_t = typename std::decay<decltype(basis_in)>::type;
        using basis_out_t = typename std::decay<decltype(basis_out)>::type;
        if constexpr (std::is_same<basis_in_t, basis_out_t>::value) {
          apply_terms(ops, basis_in, vec_in, basis_out, vec_out);
        } else {
          XDIAG_THROW(
              "Invalid combination of bases for \"ElectronDistributed\" block.")
        }
      },
      basis_in, basis_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply_dispatch(OpSum const &, ElectronDistributed const &,
                             arma::vec const &,
                             ElectronDistributed const &block, arma::vec &);
template void apply_dispatch(OpSum const &, ElectronDistributed const &,
                             arma::cx_vec const &,
                             ElectronDistributed const &block, arma::cx_vec &);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, ElectronDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    ElectronDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out) try {
  XDIAG_THROW("Apply for an OpSum on a State with a matrix not implemented yet "
              "for ElectronDistributed blocks");
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template void apply_dispatch(OpSum const &, ElectronDistributed const &,
                             arma::mat const &,
                             ElectronDistributed const &block, arma::mat &);
template void apply_dispatch(OpSum const &, ElectronDistributed const &,
                             arma::cx_mat const &,
                             ElectronDistributed const &block, arma::cx_mat &);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  using namespace basis::spinhalf_distributed;

  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      [&](auto &&basis_in, auto &&basis_out) {
        using basis_in_t = typename std::decay<decltype(basis_in)>::type;
        using basis_out_t = typename std::decay<decltype(basis_out)>::type;
        if constexpr (std::is_same<basis_in_t, basis_out_t>::value) {
          apply_terms(ops, basis_in, vec_in, basis_out, vec_out);
        } else {
          XDIAG_THROW(
              "Invalid combination of bases for \"SpinhalfDistributed\" block.")
        }
      },
      basis_in, basis_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply_dispatch(OpSum const &, SpinhalfDistributed const &,
                             arma::vec const &,
                             SpinhalfDistributed const &block, arma::vec &);
template void apply_dispatch(OpSum const &, SpinhalfDistributed const &,
                             arma::cx_vec const &,
                             SpinhalfDistributed const &block, arma::cx_vec &);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, SpinhalfDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    SpinhalfDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out) try {
  XDIAG_THROW("Apply for an OpSum on a State with a matrix not implemented yet "
              "for SpinhalfDistributed blocks");
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply_dispatch(OpSum const &, SpinhalfDistributed const &,
                             arma::mat const &,
                             SpinhalfDistributed const &block, arma::mat &);
template void apply_dispatch(OpSum const &, SpinhalfDistributed const &,
                             arma::cx_mat const &,
                             SpinhalfDistributed const &block, arma::cx_mat &);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, tJDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  using namespace basis::tj_distributed;
  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      [&](auto &&basis_in, auto &&basis_out) {
        using basis_in_t = typename std::decay<decltype(basis_in)>::type;
        using basis_out_t = typename std::decay<decltype(basis_out)>::type;
        if constexpr (std::is_same<basis_in_t, basis_out_t>::value) {
          apply_terms(ops, basis_in, vec_in, basis_out, vec_out);
        } else {
          XDIAG_THROW(
              "Invalid combination of bases for \"tJDistributed\" block.")
        }
      },
      basis_in, basis_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void apply_dispatch(OpSum const &, tJDistributed const &,
                             arma::vec const &, tJDistributed const &block,
                             arma::vec &);
template void apply_dispatch(OpSum const &, tJDistributed const &,
                             arma::cx_vec const &, tJDistributed const &block,
                             arma::cx_vec &);

template <typename coeff_t>
void apply_dispatch(OpSum const &ops, tJDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out) try {
  XDIAG_THROW("Apply for an OpSum on a State with a matrix not implemented yet "
              "for tJDistributed blocks");
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template void apply_dispatch(OpSum const &, tJDistributed const &,
                             arma::mat const &, tJDistributed const &block,
                             arma::mat &);
template void apply_dispatch(OpSum const &, tJDistributed const &,
                             arma::cx_mat const &, tJDistributed const &block,
                             arma::cx_mat &);

#endif

} // namespace xdiag::algebra
