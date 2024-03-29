#pragma once

#include <hydra/extern/armadillo/armadillo>
#include <hydra/blocks/blocks.h>
#include <hydra/common.h>

namespace hydra {

class State {
public:
  // Interface
  State() = default;

  State(block_variant_t const &block, bool real = true, int64_t n_cols = 1);

  template <typename block_t>
  explicit State(block_t const &block, bool real = true, int64_t n_cols = 1);

  template <typename block_t, typename coeff_t>
  State(block_t const &block, arma::Col<coeff_t> const &vector);

  template <typename block_t, typename coeff_t>
  State(block_t const &block, arma::Mat<coeff_t> const &matrix);

  int64_t n_sites() const;
  bool isreal() const;
  bool iscomplex() const;

  State real() const;
  State imag() const;
  void make_complex();

  int64_t dim() const;
  int64_t size() const;
  int64_t n_rows() const;
  int64_t n_cols() const;

  State col(int64_t n, bool copy = true) const;

  arma::vec vector(int64_t n = 0, bool copy = true) const;
  arma::mat matrix(bool copy = true) const;
  arma::cx_vec vectorC(int64_t n = 0, bool copy = true) const;
  arma::cx_mat matrixC(bool copy = true) const;

  // Developer section
  template <typename block_t>
  State(block_t const &block, double const *ptr, int64_t n_cols,
        int64_t stride = 1);

  template <typename block_t>
  State(block_t const &block, complex const *ptr, int64_t n_cols);
  block_variant_t block() const;

  double *memptr();
  complex *memptrC();
  double *colptr(int64_t col);
  complex *colptrC(int64_t col);
  
private:
  bool real_;

  int64_t dim_;
  int64_t n_rows_;
  int64_t n_cols_;
  block_variant_t block_;
  mutable std::vector<double> storage_;
};

State zero_state(block_variant_t const &block, bool real = true);
template <typename block_t>
State zero_state(block_t const &block, bool real = true);
State zeros_like(State const& state);
  
  
} // namespace hydra
