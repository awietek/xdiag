// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/armadillo.hpp>
#include <xdiag/blocks/blocks.hpp>

namespace xdiag {

class XDIAG_API State {
public:
  State() = default;
  explicit State(Block const &block, bool real = true, int64_t ncols = 1);
  State(Block const &block, arma::vec const &vector);
  State(Block const &block, arma::cx_vec const &vector);
  State(Block const &block, arma::mat const &matrix);
  State(Block const &block, arma::cx_mat const &matrix);

  bool isvalid() const;
  int64_t nsites() const;
  bool isreal() const;
  State real() const;
  State imag() const;
  void make_complex();
  int64_t dim() const;
  int64_t size() const;
  int64_t nrows() const;
  int64_t ncols() const;

  State col(int64_t n, bool copy = true) const;
  arma::vec vector(int64_t n = 0, bool copy = true) const;
  arma::mat matrix(bool copy = true) const;
  arma::cx_vec vectorC(int64_t n = 0, bool copy = true) const;
  arma::cx_mat matrixC(bool copy = true) const;

  // Developer section
  State(Block const &block, double const *ptr, int64_t ncols,
        int64_t stride = 1);
  State(Block const &block, complex const *ptr, int64_t ncols);
  double *memptr();
  complex *memptrC();
  double *colptr(int64_t col);
  complex *colptrC(int64_t col);
  Block block() const;

private:
  bool valid_ = false;
  Block block_;
  bool real_;
  int64_t nrows_;
  int64_t ncols_;
  mutable std::vector<double> storage_;

  void init0(bool real, int64_t nrows, int64_t ncols);
  void initcopy(const double *ptr, int64_t nrows, int64_t ncols,
                int64_t stride = 1);
  void initcopy(const complex *ptr, int64_t nrows, int64_t ncols);
};

XDIAG_API bool isvalid(State const &s);
XDIAG_API int64_t nsites(State const &s);
XDIAG_API bool isapprox(State const &v, State const &w, double rtol = 1e-12,
                        double atol = 1e-12);
XDIAG_API bool isreal(State const &s);
XDIAG_API State real(State const &s);
XDIAG_API State imag(State const &s);
XDIAG_API void make_complex(State &s);
XDIAG_API int64_t dim(State const &s);
XDIAG_API int64_t size(State const &s);
XDIAG_API int64_t nrows(State const &s);
XDIAG_API int64_t ncols(State const &s);

XDIAG_API State col(State const &s, int64_t n, bool copy = true);
XDIAG_API arma::vec vector(State const &s, int64_t n = 0, bool copy = true);
XDIAG_API arma::mat matrix(State const &s, bool copy = true);
XDIAG_API arma::cx_vec vectorC(State const &s, int64_t n = 0, bool copy = true);
XDIAG_API arma::cx_mat matrixC(State const &s, bool copy = true);

XDIAG_API std::ostream &operator<<(std::ostream &out, State const &state);
XDIAG_API std::string to_string(State const &state);

// arithmetic operators
XDIAG_API State &operator*=(State &X, double alpha);
XDIAG_API State operator*(State const &X, double alpha);
XDIAG_API State operator*(double alpha, State const &X);

XDIAG_API State &operator*=(State &X, complex alpha);
XDIAG_API State operator*(State const &X, complex alpha);
XDIAG_API State operator*(complex alpha, State const &X);

XDIAG_API State &operator/=(State &X, double alpha);
XDIAG_API State operator/(State const &X, double alpha);

XDIAG_API State &operator/=(State &X, complex alpha);
XDIAG_API State operator/(State const &X, complex alpha);

XDIAG_API State &operator+=(State &v, State const &w);
XDIAG_API State &operator-=(State &v, State const &w);
XDIAG_API State operator+(State const &v, State const &w);
XDIAG_API State operator-(State const &v, State const &w);

XDIAG_API State operator-(State const &v);
} // namespace xdiag
