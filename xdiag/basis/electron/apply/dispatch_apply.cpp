#include "dispatch_apply.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/basis/electron/apply/dispatch.hpp>

namespace xdiag::basis::electron {

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Electron const &block_in,
                    arma::Col<coeff_t> const &vec_in, Electron const &block_out,
                    arma::Col<coeff_t> &vec_out) try {
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(vec_in, vec_out, idx_in, idx_out, val);
  };
  dispatch<coeff_t>(ops, block_in, block_out, fill);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void dispatch_apply(OpSum const &, Electron const &, arma::vec const &,
                             Electron const &block, arma::vec &);
template void dispatch_apply(OpSum const &, Electron const &,
                             arma::cx_vec const &, Electron const &block,
                             arma::cx_vec &);

template <typename coeff_t>
void dispatch_apply(OpSum const &ops, Electron const &block_in,
                    arma::Mat<coeff_t> const &mat_in, Electron const &block_out,
                    arma::Mat<coeff_t> &mat_out) try {
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    return fill_apply(mat_in, mat_out, idx_in, idx_out, val);
  };
  dispatch<coeff_t>(ops, block_in, block_out, fill);
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

template void dispatch_apply(OpSum const &, Electron const &, arma::mat const &,
                             Electron const &block, arma::mat &);

template void dispatch_apply(OpSum const &, Electron const &, arma::cx_mat const &,
                             Electron const &block, arma::cx_mat &);

} // namespace xdiag::basis::electron
