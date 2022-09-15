#include "gpwf_spinhalf.h"

#include <hydra/bitops/bitops.h>

namespace hydra {

template <class coeff_t>
GPWFSpinhalf<coeff_t>::GPWFSpinhalf(int n_sites, arma::Mat<coeff_t> onebody_wfs,
                                    int n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_),
      work_matrix_(arma::zeros<arma::Mat<coeff_t>>(n_sites, n_sites)) {
  using arma::span;
  onebody_wfs_up_ = onebody_wfs(span(0, n_sites), span(0, n_up_));
  onebody_wfs_dn_ = onebody_wfs(span(0, n_sites), span(n_up_, n_sites_));
  assert((int)onebody_wfs.n_rows == (int)n_sites);
  assert((int)onebody_wfs.n_cols == (int)n_sites);
}

template <class coeff_t>
coeff_t GPWFSpinhalf<coeff_t>::coefficient(uint64_t state, bool print_work) {

  work_matrix_.zeros();
  for (int i = 0; i < n_sites_; ++i) {

    // spin at site i is up
    if (bitops::gbit(state, i)) {
      // for (int k = 0; k < n_up_; ++k)
      //   work_matrix_(k, i) = onebody_wfs_up_(i, k);
      work_matrix_(arma::span(0, n_up_), i) =
          onebody_wfs_up_(i, arma::span(0, n_up_));
    }

    // spin at site i is dn
    else {
      // for (int k = 0; k < n_dn_; ++k)
      //   work_matrix_(k + n_up_, i) = onebody_wfs_dn_(i, k);
      work_matrix_(arma::span(n_up_, n_sites_), i) =
          onebody_wfs_dn_(i, arma::span(0, n_dn_));
    }
  }
  // if (print_work) {
  //   std::cout << n_sites_ << " " << n_up_ << " " << n_dn_ << "\n";
  //   LilaPrint(Real(onebody_wfs_up_));
  //   LilaPrint(Real(onebody_wfs_dn_));
  //   LilaPrint(Real(work_matrix_));
  //   exit(1);
  // }
  return arma::det(work_matrix_);
}

template <class coeff_t>
bool GPWFSpinhalf<coeff_t>::operator==(GPWFSpinhalf<coeff_t> const &other) {
  return ((n_sites_ == other.n_sites_) && (n_up_ == other.n_up_) &&
          (n_dn_ == other.n_dn_) &&
          (arma::approx_equal(onebody_wfs_up_, other.onebody_wfs_up_, "both", 1e-12, 1e-12)) &&
          (arma::approx_equal(onebody_wfs_dn_, other.onebody_wfs_dn_, "both", 1e-12, 1e-12)));
}

template class GPWFSpinhalf<double>;
template class GPWFSpinhalf<complex>;

} // namespace hydra
