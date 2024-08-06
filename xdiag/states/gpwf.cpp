#include "gpwf.hpp"

#include <xdiag/bits/bitops.hpp>

namespace xdiag {

GPWF::GPWF(int64_t n_sites, arma::mat onebody_wfs, int64_t n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_), isreal_(true) {
  using namespace arma;
  work_matrix_ = mat(n_sites, n_sites, fill::zeros);
  onebody_wfs_up_ = onebody_wfs(span(0, n_sites), span(0, n_up_));
  onebody_wfs_dn_ = onebody_wfs(span(0, n_sites), span(n_up_, n_sites));

  work_matrix_c_ = cx_mat();
  onebody_wfs_up_c_ = cx_mat();
  onebody_wfs_dn_c_ = cx_mat();
}

GPWF::GPWF(int64_t n_sites, arma::cx_mat onebody_wfs, int64_t n_up)
    : n_sites_(n_sites), n_up_(n_up >= 0 ? n_up : n_sites / 2),
      n_dn_(n_sites - n_up_), isreal_(false) {
  using namespace arma;
  work_matrix_ = mat();
  onebody_wfs_up_ = mat();
  onebody_wfs_dn_ = mat();

  work_matrix_c_ = cx_mat(n_sites, n_sites, fill::zeros);
  onebody_wfs_up_c_ = onebody_wfs(span(0, n_sites), span(0, n_up_));
  onebody_wfs_dn_c_ = onebody_wfs(span(0, n_sites), span(n_up_, n_sites));
}

int64_t GPWF::n_sites() const { return n_sites_; }
int64_t GPWF::n_up() const { return n_up_; }
bool GPWF::isreal() const { return isreal_; }

template <class coeff_t>
static coeff_t gpwf_coefficient(ProductState const &pstate, int64_t n_sites,
                                int64_t n_up, int64_t n_dn,
                                arma::Mat<coeff_t> const &onebody_wfs_up,
                                arma::Mat<coeff_t> const &onebody_wfs_dn,
                                arma::Mat<coeff_t> &work_matrix) {
  work_matrix.zeros();
  for (int64_t i = 0; i < n_sites; ++i) {
    std::string s = pstate[i];
    if (s == "Up") {
      work_matrix(arma::span(0, n_up), i) =
          onebody_wfs_up(i, arma::span(0, n_up));
    } else { // spin at site i is dn
      work_matrix(arma::span(n_up, n_sites), i) =
          onebody_wfs_dn(i, arma::span(0, n_dn));
    }
  }
  return arma::det(work_matrix);
}

double GPWF::coefficient(ProductState const &pstate) const try {
  if (pstate.n_sites() != n_sites_) {
    XDIAG_THROW("Number of sites different between ProductState and GPWF");
  }
  if (!isreal_) {
    XDIAG_THROW(
        "Cannot compute a real coefficient of a genuinely complex GPWF");
  }
  return gpwf_coefficient(pstate, n_sites_, n_up_, n_dn_, onebody_wfs_up_,
                          onebody_wfs_dn_, work_matrix_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

complex GPWF::coefficientC(ProductState const &pstate) const try {
  if (pstate.n_sites() != n_sites_) {
    XDIAG_THROW("Number of sites different between ProductState and GPWF");
  }
  if (isreal_) {
    return (complex)gpwf_coefficient(pstate, n_sites_, n_up_, n_dn_,
                                     onebody_wfs_up_, onebody_wfs_dn_,
                                     work_matrix_);
  } else {
    return gpwf_coefficient(pstate, n_sites_, n_up_, n_dn_, onebody_wfs_up_c_,
                            onebody_wfs_dn_c_, work_matrix_c_);
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool GPWF::operator==(GPWF const &rhs) const {
  if ((n_sites_ != rhs.n_sites_) || (n_up_ != rhs.n_up_) ||
      (n_dn_ != rhs.n_dn_) || (isreal_ != rhs.isreal_)) {
    return false;
  } else {
    if (isreal_) {
      return (arma::approx_equal(onebody_wfs_up_, rhs.onebody_wfs_up_, "both",
                                 1e-12, 1e-12)) &&
             (arma::approx_equal(onebody_wfs_dn_, rhs.onebody_wfs_dn_, "both",
                                 1e-12, 1e-12));
    } else {
      return (arma::approx_equal(onebody_wfs_up_c_, rhs.onebody_wfs_up_c_,
                                 "both", 1e-12, 1e-12)) &&
             (arma::approx_equal(onebody_wfs_dn_c_, rhs.onebody_wfs_dn_c_,
                                 "both", 1e-12, 1e-12));
    }
  }
}

bool GPWF::operator!=(GPWF const &rhs) const { return !operator==(rhs); }

} // namespace xdiag
