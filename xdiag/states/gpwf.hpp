// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>
#include <xdiag/states/product_state.hpp>

namespace xdiag {

class GPWF {
public:
  GPWF() = default;
  GPWF(arma::mat const &onebody_wfs, int64_t nup = -1);
  GPWF(arma::cx_mat const &onebody_wfs, int64_t nup = -1);

  int64_t nsites() const;
  int64_t nup() const;
  bool isreal() const;

  double coefficient(ProductState const &state) const;
  complex coefficientC(ProductState const &state) const;

  bool operator==(GPWF const &other) const;
  bool operator!=(GPWF const &other) const;

private:
  bool isreal_;
  int64_t nsites_;
  int64_t nup_;
  int64_t ndn_;

  // Real version
  mutable arma::mat work_matrix_;
  arma::mat onebody_wfs_up_;
  arma::mat onebody_wfs_dn_;

  // Complex version
  mutable arma::cx_mat work_matrix_c_;
  arma::cx_mat onebody_wfs_up_c_;
  arma::cx_mat onebody_wfs_dn_c_;

  friend std::ostream &operator<<(std::ostream &os, GPWF const &p);
};

std::ostream &operator<<(std::ostream &out, GPWF const &state);
std::string to_string(GPWF const &state);

} // namespace xdiag
