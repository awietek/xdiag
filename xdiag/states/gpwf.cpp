#include "gpwf.hpp"

#include <xdiag/bits/bitops.hpp>

namespace xdiag {

GPWF::GPWF(arma::mat const &onebody_wfs, int64_t nup) try : isreal_(true) {
  using namespace arma;

  if (onebody_wfs.n_rows != onebody_wfs.n_cols) {
    XDIAG_THROW("Matrix for onebody wave functions must be square");
  }

  nsites_ = onebody_wfs.n_rows;
  nup_ = nup >= 0 ? nup : nsites_ / 2;
  ndn_ = nsites_ - nup_;

  work_matrix_ = mat(nsites_, nsites_, fill::zeros);
  onebody_wfs_up_ = onebody_wfs(span::all, span(0, nup_-1));
  onebody_wfs_dn_ = onebody_wfs(span::all, span(nup_, nsites_-1));

  work_matrix_c_ = cx_mat();
  onebody_wfs_up_c_ = cx_mat();
  onebody_wfs_dn_c_ = cx_mat();
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

GPWF::GPWF(arma::cx_mat const &onebody_wfs, int64_t nup) try : isreal_(false) {
  using namespace arma;

  if (onebody_wfs.n_rows != onebody_wfs.n_cols) {
    XDIAG_THROW("Matrix for onebody wave functions must be square");
  }

  nsites_ = onebody_wfs.n_rows;
  nup_ = nup >= 0 ? nup : nsites_ / 2;
  ndn_ = nsites_ - nup_;

  work_matrix_ = mat();
  onebody_wfs_up_ = mat();
  onebody_wfs_dn_ = mat();

  work_matrix_c_ = cx_mat(nsites_, nsites_, fill::zeros);
  onebody_wfs_up_c_ = onebody_wfs(span::all, span(0, nup_-1));
  onebody_wfs_dn_c_ = onebody_wfs(span::all, span(nup_, nsites_-1));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t GPWF::nsites() const { return nsites_; }
int64_t GPWF::nup() const { return nup_; }
bool GPWF::isreal() const { return isreal_; }

template <class coeff_t>
static coeff_t gpwf_coefficient(ProductState const &pstate, int64_t nsites,
                                int64_t nup, int64_t ndn,
                                arma::Mat<coeff_t> const &onebody_wfs_up,
                                arma::Mat<coeff_t> const &onebody_wfs_dn,
                                arma::Mat<coeff_t> &work_matrix) {
  work_matrix.zeros();
  for (int64_t i = 0; i < nsites; ++i) {
    std::string s = pstate[i];
    if (s == "Up") {
      work_matrix(i, arma::span(0, nup-1)) =
          onebody_wfs_up(i, arma::span(0, nup-1));
    } else { // spin at site i is dn
      work_matrix(i, arma::span(nup, nsites-1)) =
          onebody_wfs_dn(i, arma::span(0, ndn-1));
    }
  }
  return arma::det(work_matrix);
}

double GPWF::coefficient(ProductState const &pstate) const try {
  if (pstate.nsites() != nsites_) {
    XDIAG_THROW("Number of sites different between ProductState and GPWF");
  }
  if (!isreal_) {
    XDIAG_THROW(
        "Cannot compute a real coefficient of a genuinely complex GPWF");
  }
  return gpwf_coefficient(pstate, nsites_, nup_, ndn_, onebody_wfs_up_,
                          onebody_wfs_dn_, work_matrix_);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

complex GPWF::coefficientC(ProductState const &pstate) const try {
  if (pstate.nsites() != nsites_) {
    XDIAG_THROW("Number of sites different between ProductState and GPWF");
  }
  if (isreal_) {
    return (complex)gpwf_coefficient(pstate, nsites_, nup_, ndn_,
                                     onebody_wfs_up_, onebody_wfs_dn_,
                                     work_matrix_);
  } else {
    return gpwf_coefficient(pstate, nsites_, nup_, ndn_, onebody_wfs_up_c_,
                            onebody_wfs_dn_c_, work_matrix_c_);
  }

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool GPWF::operator==(GPWF const &rhs) const {
  if ((nsites_ != rhs.nsites_) || (nup_ != rhs.nup_) ||
      (ndn_ != rhs.ndn_) || (isreal_ != rhs.isreal_)) {
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

std::ostream &operator<<(std::ostream &out, GPWF const &state) {
  if (state.isreal()) {
    out << "Gutzwiller Projected WF: \n"
        << " ups:\n";
    state.onebody_wfs_up_.brief_print(out);
    out << " dns:\n";
    state.onebody_wfs_dn_.brief_print(out);
  } else {
    out << "Gutzwiller Projected WF: \n"
        << " ups:\n";
    state.onebody_wfs_up_c_.brief_print(out);
    out << " dns:\n";
    state.onebody_wfs_dn_c_.brief_print(out);
  }
  return out;
}
std::string to_string(GPWF const &state) { return to_string_generic(state); }

} // namespace xdiag
