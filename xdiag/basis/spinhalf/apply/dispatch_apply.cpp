#include "dispatch_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch.hpp>

namespace xdiag::basis {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Spinhalf const &block_in,
                    arma::Col<coeff_t> const &vec_in, Spinhalf const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  spinhalf::dispatch<coeff_t>(ops, block_in, block_out, fill);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Spinhalf const &block_in,
                    arma::Mat<coeff_t> const &mat_in, Spinhalf const &block_out,
                    arma::Mat<coeff_t> &mat_out) try {
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(mat_in, mat_out, idx_in, idx_out, val);
  };
  spinhalf::dispatch<coeff_t>(ops, block_in, block_out, fill);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void dispatch_apply(OpSum const &, Spinhalf const &, arma::vec const &,
                             Spinhalf const &block, arma::vec &);
template void dispatch_apply(OpSum const &, Spinhalf const &,
                             arma::cx_vec const &, Spinhalf const &block,
                             arma::cx_vec &);
template void dispatch_apply(OpSum const &, Spinhalf const &, arma::mat const &,
                             Spinhalf const &block, arma::mat &);
template void dispatch_apply(OpSum const &, Spinhalf const &,
                             arma::cx_mat const &, Spinhalf const &block,
                             arma::cx_mat &);

} // namespace xdiag::basis
