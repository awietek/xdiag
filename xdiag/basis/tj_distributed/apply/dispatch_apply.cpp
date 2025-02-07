#include "dispatch_apply.hpp"

#include <xdiag/basis/tj_distributed/apply/apply_terms.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJDistributed const &block_in,
                    arma::Col<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  auto const &basis_in = block_in.basis();
  auto const &basis_out = block_out.basis();

  std::visit(
      [&](auto &&basis_in, auto &&basis_out) {
        using basis_in_t = typename std::decay<decltype(basis_in)>::type;
        using basis_out_t = typename std::decay<decltype(basis_out)>::type;
        if constexpr (std::is_same<basis_in_t, basis_out_t>::value) {
          basis::tj_distributed::apply_terms(ops, basis_in, vec_in, basis_out,
                                             vec_out);
        } else {
          XDIAG_THROW(
              "Invalid combination of bases for \"tJDistributed\" block.")
        }
      },
      basis_in, basis_out);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void dispatch_apply(OpSum const &, tJDistributed const &,
                             arma::vec const &, tJDistributed const &block,
                             arma::vec &);
template void dispatch_apply(OpSum const &, tJDistributed const &,
                             arma::cx_vec const &, tJDistributed const &block,
                             arma::cx_vec &);

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, tJDistributed const &block_in,
                    arma::Mat<coeff_t> const &vec_in,
                    tJDistributed const &block_out,
                    arma::Mat<coeff_t> &vec_out) try {
  XDIAG_THROW("Apply for an OpSum on a State with a matrix not implemented yet "
              "for tJDistributed blocks");
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
template void dispatch_apply(OpSum const &, tJDistributed const &,
                             arma::mat const &, tJDistributed const &block,
                             arma::mat &);
template void dispatch_apply(OpSum const &, tJDistributed const &,
                             arma::cx_mat const &, tJDistributed const &block,
                             arma::cx_mat &);

} // namespace xdiag::basis
