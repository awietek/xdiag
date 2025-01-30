#pragma once

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/common.hpp>
#include <xdiag/extern/armadillo/armadillo>

namespace xdiag {

class State {
public:
  XDIAG_API State() = default;
  XDIAG_API explicit State(Block const &block, bool real = true,
                           int64_t ncols = 1);

  template <typename block_t>
  XDIAG_API explicit State(block_t const &block, bool real = true,
                           int64_t ncols = 1);

  template <typename block_t, typename coeff_t>
  XDIAG_API State(block_t const &block, arma::Col<coeff_t> const &vector);

  template <typename block_t, typename coeff_t>
  XDIAG_API State(block_t const &block, arma::Mat<coeff_t> const &matrix);

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
  template <typename block_t>
  State(block_t const &block, double const *ptr, int64_t ncols,
        int64_t stride = 1);

  template <typename block_t>
  State(block_t const &block, complex const *ptr, int64_t ncols);
  Block block() const;

  double *memptr();
  complex *memptrC();
  double *colptr(int64_t col);
  complex *colptrC(int64_t col);

private:
  bool real_;
  int64_t dim_;
  int64_t nrows_;
  int64_t ncols_;
  Block block_;
  mutable std::vector<double> storage_;
};

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

} // namespace xdiag
