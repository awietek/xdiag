#include "state.h"

#include <xdiag/blocks/blocks.h>

namespace xdiag {

State::State(block_variant_t const &block, bool real, int64_t n_cols)
    : real_(real), dim_(xdiag::dim(block)), n_rows_(xdiag::size(block)),
      n_cols_(n_cols), block_(block) {
  try {
    if (real) {
      resize_vector(storage_, size());
    } else {
      resize_vector(storage_, 2 * size());
    }
  } catch (...) {
    XDiagRethrow("Unable to allocate memory for State");
  }
}

template <typename block_t>
State::State(block_t const &block, bool real, int64_t n_cols)
    : real_(real), dim_(block.dim()), n_rows_(block.size()), n_cols_(n_cols),
      block_(block) {
  try {
    if (real) {
      resize_vector(storage_, size());
    } else {
      resize_vector(storage_, 2 * size());
    }
  } catch (...) {
    XDiagRethrow("Unable to allocate memory for State");
  }
}

template <typename block_t>
State::State(block_t const &block, double const *ptr, int64_t n_cols,
             int64_t stride)
    : real_(true), n_rows_(block.size()), n_cols_(n_cols), block_(block) {

  try {
    resize_vector(storage_, size());
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to allocate memory for State");
  }

  try {
    if (stride == 1) {
      std::copy(ptr, ptr + size(), storage_.data());
    } else {
      for (int64_t i = 0, is = 0; i < size(); ++i, is += stride) {
        storage_[i] = ptr[is];
      }
    }
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to copy memory for State");
  }
}

template <typename block_t>
State::State(block_t const &block, complex const *ptr, int64_t n_cols)
    : real_(false), n_rows_(block.size()), n_cols_(n_cols), block_(block) {
  try {
    resize_vector(storage_, 2 * size());
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to allocate memory for State");
  }

  try {
    std::copy(ptr, ptr + size(), reinterpret_cast<complex *>(storage_.data()));
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to copy memory for State");
  }
}

template <typename block_t, typename coeff_t>
State::State(block_t const &block, arma::Col<coeff_t> const &vector)
    : real_(xdiag::isreal<coeff_t>()), n_rows_(block.size()), n_cols_(1),
      block_(block) {
  if (block.size() != (int64_t)vector.n_rows) {
    XDiagThrow(std::invalid_argument,
               "Block dimension not equal to vector dimension");
  }

  try {
    if (real_) {
      resize_vector(storage_, size());
    } else {
      resize_vector(storage_, 2 * size());
    }
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to allocate memory for State");
  }

  try {
    std::copy(vector.memptr(), vector.memptr() + size(),
              reinterpret_cast<coeff_t *>(storage_.data()));
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to copy memory for State");
  }
}

template <typename block_t, typename coeff_t>
State::State(block_t const &block, arma::Mat<coeff_t> const &matrix)
    : real_(xdiag::isreal<coeff_t>()), n_rows_(matrix.n_rows),
      n_cols_(matrix.n_cols), block_(block) {

  if (block.size() != (int64_t)matrix.n_rows) {
    XDiagThrow(std::invalid_argument,
               "Block dimension not equal to number of rows in matrix");
  }

  try {
    if (real_) {
      resize_vector(storage_, size());
    } else {
      resize_vector(storage_, 2 * size());
    }
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to allocate memory for State");
  }

  try {
    std::copy(matrix.memptr(), matrix.memptr() + size(),
              reinterpret_cast<coeff_t *>(storage_.data()));
  } catch (...) {
    XDiagThrow(std::runtime_error, "Unable to copy memory for State");
  }
}

int64_t State::n_sites() const { return xdiag::n_sites(block_); }
bool State::isreal() const { return real_; }
bool State::iscomplex() const { return !isreal(); }

State State::real() const {
  if (isreal()) {
    return (*this);
  } else {
    // Return a State with zero values
    double *ptr = storage_.data();
    return std::visit(
        [&](auto &&block) { return State(block, ptr, n_cols_, 2); }, block_);
  }
}

State State::imag() const {
  if (isreal()) {
    return (*this);
  } else {
    return State(block_, true, n_cols_);
  }
}

void State::make_complex() try {
  if (isreal()) {
    real_ = false;

    try {
      resize_vector(storage_, 2 * size());
    } catch (...) {
      XDiagThrow(std::runtime_error, "Unable to allocate memory for State");
    }

    double *ptr = storage_.data();
    for (int64_t i = size() - 1; i >= 0; --i) {
      std::swap(ptr[i << 1], ptr[i]);
    }
  }
} catch (...) {
  XDiagRethrow("Error transforming to complex State");
}

int64_t State::dim() const { return xdiag::dim(block_); }

int64_t State::size() const { return n_rows_ * n_cols_; }
int64_t State::n_rows() const { return n_rows_; }
int64_t State::n_cols() const { return n_cols_; }
block_variant_t State::block() const { return block_; }

State State::col(int64_t n, bool copy) const try {
  if (n >= n_cols_) {
    XDiagThrow(std::invalid_argument,
               "Column index larger than the number of columns");
  }
  if (real_) {
    return std::visit([&](auto &&b) { return State(b, vector(n, copy)); },
                      block_);
  } else {
    return std::visit([&](auto &&b) { return State(b, vectorC(n, copy)); },
                      block_);
  }
} catch (...) {
  XDiagRethrow("Unable to retrieve column State");
  return State();
}

arma::vec State::vector(int64_t n, bool copy) const try {
  if (!real_) {
    XDiagThrow(std::logic_error,
               "Cannot return a real armadillo vector from a "
               "complex State (maybe use vectorC(...) instead)");
  } else if (n >= n_cols_) {
    XDiagThrow(std::invalid_argument,
               "Column index larger than the number of columns");
  } else if (n < 0) {
    XDiagThrow(std::invalid_argument, "Negative column index");
  }
  return arma::vec(storage_.data() + n * n_rows_, n_rows_, copy, !copy);
} catch (...) {
  XDiagRethrow("Unable to access real vector of State");
  return arma::vec();
}

arma::mat State::matrix(bool copy) const try {
  if (!real_) {
    XDiagThrow(std::logic_error,
               "Cannot return a real armadillo matrix from a "
               "complex state (maybe use matrixC(...) instead)");
  }
  return arma::mat(storage_.data(), n_rows_, n_cols_, copy, !copy);
} catch (...) {
  XDiagRethrow("Unable to access real matrix of State");
  return arma::mat();
}

arma::cx_vec State::vectorC(int64_t n, bool copy) const try {
  if (real_) {
    XDiagThrow(std::logic_error,
               "Cannot return a complex armadillo vector from a "
               "real State (maybe use vector(...) instead)");
  } else if (n >= n_cols_) {
    XDiagThrow(std::invalid_argument,
               "Column index larger than the number of columns");
  } else if (n < 0) {
    XDiagThrow(std::invalid_argument, "Negative column index");
  }
  return arma::cx_vec(reinterpret_cast<complex *>(storage_.data()) +
                          n * n_rows_,
                      n_rows_, copy, !copy);
} catch (...) {
  XDiagRethrow("Unable to access complex vector of State");
  return arma::cx_vec();
}

arma::cx_mat State::matrixC(bool copy) const try {
  if (real_) {
    XDiagThrow(std::logic_error,
               "Cannot return a complex armadillo matrix from a "
               "real state (maybe use matrix(...) instead)");
  }
  return arma::cx_mat(reinterpret_cast<complex *>(storage_.data()), n_rows_,
                      n_cols_, copy, !copy);
} catch (...) {
  XDiagRethrow("Unable to access complex matrix of State");
  return arma::cx_mat();
}

double *State::memptr() { return storage_.data(); }
complex *State::memptrC() {
  return reinterpret_cast<complex *>(storage_.data());
}
double *State::colptr(int64_t col) {
  if ((col < 0) || (col >= n_cols_)) {
    XDiagThrow(std::invalid_argument, "Invalid column index requested");
  }
  return memptr() + col * n_rows_;
}
complex *State::colptrC(int64_t col) {
  if ((col < 0) || (col >= n_cols_)) {
    XDiagThrow(std::invalid_argument, "Invalid column index requested");
  }
  return memptrC() + col * n_rows_;
}

State zeros_like(State const &state) {
  return State(state.block(), state.isreal(), state.n_cols());
}

template State::State(Spinhalf const &, bool, int64_t);
template State::State(tJ const &, bool, int64_t);
template State::State(Electron const &, bool, int64_t);

template State::State(Spinhalf const &, double const *, int64_t, int64_t);
template State::State(tJ const &, double const *, int64_t, int64_t);
template State::State(Electron const &, double const *, int64_t, int64_t);

template State::State(Spinhalf const &block, complex const *, int64_t size);
template State::State(tJ const &block, complex const *, int64_t size);
template State::State(Electron const &block, complex const *, int64_t size);

template State::State(Spinhalf const &block, arma::Col<double> const &vector);
template State::State(tJ const &block, arma::Col<double> const &vector);
template State::State(Electron const &block, arma::Col<double> const &vector);
template State::State(Spinhalf const &block, arma::Col<complex> const &vector);
template State::State(tJ const &block, arma::Col<complex> const &vector);
template State::State(Electron const &block, arma::Col<complex> const &vector);

template State::State(Spinhalf const &block, arma::Mat<double> const &vector);
template State::State(tJ const &block, arma::Mat<double> const &vector);
template State::State(Electron const &block, arma::Mat<double> const &vector);
template State::State(Spinhalf const &block, arma::Mat<complex> const &vector);
template State::State(tJ const &block, arma::Mat<complex> const &vector);
template State::State(Electron const &block, arma::Mat<complex> const &vector);

#ifdef XDIAG_USE_MPI
template State::State(tJDistributed const &, bool, int64_t);
template State::State(tJDistributed const &, double const *, int64_t, int64_t);
template State::State(tJDistributed const &block, complex const *ptr,
                      int64_t size);
template State::State(tJDistributed const &block,
                      arma::Col<double> const &vector);
template State::State(tJDistributed const &block,
                      arma::Col<complex> const &vector);
template State::State(tJDistributed const &block,
                      arma::Mat<double> const &vector);
template State::State(tJDistributed const &block,
                      arma::Mat<complex> const &vector);

#endif

} // namespace xdiag
