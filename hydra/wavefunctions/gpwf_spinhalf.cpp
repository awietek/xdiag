#include "gpwf_spinhalf.h"

#include <hydra/utils/bitops.h>

namespace hydra {

using namespace lila;

template <class bit_t, class coeff_t>
GPWFSpinhalf<bit_t, coeff_t>::GPWFSpinhalf(
    int n_sites, lila::Matrix<coeff_t> const &onebody_wfs, int n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_), work_matrix_(Zeros(2 * n_sites, 2 * n_sites)) {
  assert(onebody_wfs.m() == n_sites);
  assert(onebody_wfs.n() == n_sites);

  onebody_wfs_up_ = onebody_wfs(ALL, {0, n_up_});
  onebody_wfs_dn_ = onebody_wfs(ALL, {n_up_, n_sites});
}

coeff_t template <class bit_t, class coeff_t>
GPWFSpinhalf<bit_t, coeff_t>::coefficient(bit_t const &state) {

  Zeros(work_matrix_);
  for (int i = 0; i < n_sites; ++i) {

    // spin at site i is up
    if (gbit(state, i)) {
      // for (int k = 0; k < n_ups; ++k)
      //   work_matrix_(k, i) = onebody_wfs_up_(k, i);
      work_matrix_({0 : n_up_}, i) = onebody_wfs_up_(i, {0, n_up_});
    }

    // spin at site i is dn
    else {
      // for (int k = 0; k < n_dns_; ++k)
      //   work_matrix_(k + n_up_, i) = onebody_wfs_dn_(k, i);
      work_matrix_({n_up_, n_sites}, i) = onebody_wfs_dn(i, {0, n_dn_})
    }
  }
  return DeterminantInplace(work_matrix_);
}

} // namespace hydra
