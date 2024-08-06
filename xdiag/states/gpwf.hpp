#pragma once

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class GPWF {
public:
  GPWF() = default;
  GPWF(int64_t n_sites, arma::mat onebody_wfs, int64_t n_up = -1);
  GPWF(int64_t n_sites, arma::cx_mat onebody_wfs, int64_t n_up = -1);

  int64_t n_sites() const;
  int64_t n_up() const;
  bool isreal() const;

  double coefficient(ProductState const &state) const;
  complex coefficientC(ProductState const &state) const;

  bool operator==(GPWF const &other) const;
  bool operator!=(GPWF const &other) const;

private:
  int64_t n_sites_;
  int64_t n_up_;
  int64_t n_dn_;
  bool isreal_;

  // Real version
  mutable arma::mat work_matrix_;
  arma::mat onebody_wfs_up_;
  arma::mat onebody_wfs_dn_;

  // Complex version
  mutable arma::cx_mat work_matrix_c_;
  arma::cx_mat onebody_wfs_up_c_;
  arma::cx_mat onebody_wfs_dn_c_;
};

} // namespace xdiag
