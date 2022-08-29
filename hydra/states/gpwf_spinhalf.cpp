#include "gpwf_spinhalf.h"

#include <hydra/bitops/bitops.h>
#include <hydra/utils/complex.h>

namespace hydra {

using namespace lila;

template <class coeff_t>
GPWFSpinhalf<coeff_t>::GPWFSpinhalf(int n_sites,
                                    lila::Matrix<coeff_t> onebody_wfs, int n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_),
      work_matrix_(lila::Zeros<coeff_t>(n_sites, n_sites)) {
  onebody_wfs_up_ = onebody_wfs({0, n_sites}, {0, n_up_});
  onebody_wfs_dn_ = onebody_wfs({0, n_sites}, {n_up_, n_sites_});
  assert(onebody_wfs.m() == n_sites);
  assert(onebody_wfs.n() == n_sites);
}

template <class coeff_t>
coeff_t GPWFSpinhalf<coeff_t>::coefficient(uint64_t state, bool print_work) {

  lila::Zeros(work_matrix_);
  for (int i = 0; i < n_sites_; ++i) {

    // spin at site i is up
    if (bitops::gbit(state, i)) {
      // for (int k = 0; k < n_up_; ++k)
      //   work_matrix_(k, i) = onebody_wfs_up_(i, k);
      work_matrix_({0, n_up_}, i) = onebody_wfs_up_(i, {0, n_up_});
    }

    // spin at site i is dn
    else {
      // for (int k = 0; k < n_dn_; ++k)
      //   work_matrix_(k + n_up_, i) = onebody_wfs_dn_(i, k);
      work_matrix_({n_up_, n_sites_}, i) = onebody_wfs_dn_(i, {0, n_dn_});
    }
  }
  if (print_work) {
    std::cout << n_sites_ << " " << n_up_ << " " << n_dn_ << "\n";
    LilaPrint(Real(onebody_wfs_up_));
    LilaPrint(Real(onebody_wfs_dn_));
    LilaPrint(Real(work_matrix_));
    exit(1);
  }
  return DeterminantInplace(work_matrix_);
}

template <class coeff_t>
bool GPWFSpinhalf<coeff_t>::operator==(GPWFSpinhalf<coeff_t> const &other) {
  return ((n_sites_ == other.n_sites_) && (n_up_ == other.n_up_) &&
          (n_dn_ == other.n_dn_) &&
          (onebody_wfs_up_ == other.onebody_wfs_up_) &&
          (onebody_wfs_dn_ == other.onebody_wfs_dn_));
}

template class GPWFSpinhalf<double>;
template class GPWFSpinhalf<complex>;

} // namespace hydra
