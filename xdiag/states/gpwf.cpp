#include "gpwf.h"
#include <xdiag/bits/bitops.h>

namespace xdiag {

template <class coeff_t>
GPWF<coeff_t>::GPWF(int64_t n_sites, arma::Mat<coeff_t> onebody_wfs,
                    int64_t n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_),
      work_matrix_(arma::zeros<arma::Mat<coeff_t>>(n_sites, n_sites)) {
  using arma::span;
  onebody_wfs_up_ = onebody_wfs(span(0, n_sites), span(0, n_up_));
  onebody_wfs_dn_ = onebody_wfs(span(0, n_sites), span(n_up_, n_sites_));
  assert((int64_t)onebody_wfs.n_rows == (int64_t)n_sites);
  assert((int64_t)onebody_wfs.n_cols == (int64_t)n_sites);
}

template <class coeff_t> coeff_t GPWF<coeff_t>::coefficient(uint64_t state) {
  work_matrix_.zeros();
  for (int64_t i = 0; i < n_sites_; ++i) {
    // spin at site i is up
    if (bits::gbit(state, i)) {
      work_matrix_(arma::span(0, n_up_), i) =
          onebody_wfs_up_(i, arma::span(0, n_up_));
    }

    // spin at site i is dn
    else {
      work_matrix_(arma::span(n_up_, n_sites_), i) =
          onebody_wfs_dn_(i, arma::span(0, n_dn_));
    }
  }
  return arma::det(work_matrix_);
}

template <class coeff_t>
bool GPWF<coeff_t>::operator==(GPWF<coeff_t> const &other) {
  return ((n_sites_ == other.n_sites_) && (n_up_ == other.n_up_) &&
          (n_dn_ == other.n_dn_) &&
          (arma::approx_equal(onebody_wfs_up_, other.onebody_wfs_up_, "both",
                              1e-12, 1e-12)) &&
          (arma::approx_equal(onebody_wfs_dn_, other.onebody_wfs_dn_, "both",
                              1e-12, 1e-12)));
}

template class GPWF<double>;
template class GPWF<complex>;

} // namespace xdiag
