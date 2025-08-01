// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/utils/arma_to_cx.hpp>

namespace xdiag {

static void safe_resize(std::vector<double> &vec, int64_t size) try {
  vec.resize(size, 0);
} catch (...) {
  XDIAG_THROW("Unable to resize vector");
}

// This initialization just allocates the memory
void State::init0(bool real, int64_t nrows, int64_t ncols) try {
  if (ncols < 1) {
    XDIAG_THROW(fmt::format("Invalid number of columns for State given: {}. "
                            "Needs to be an integer larger equal to 1.",
                            ncols));
  }
  real_ = real;
  nrows_ = nrows;
  ncols_ = ncols;
  if (real) {
    safe_resize(storage_, size());
  } else {
    safe_resize(storage_, 2 * size());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// This initialization copies the memory (double)
void State::initcopy(const double *ptr, int64_t nrows, int64_t ncols,
                     int64_t stride) try {
  init0(true, nrows, ncols);
  if (stride == 1) {
    std::copy(ptr, ptr + size(), storage_.data());
  } else if (stride == 2) {
    for (int64_t i = 0, is = 0; i < size(); ++i, is += 2) {
      storage_[i] = ptr[is];
    }
  } else {
    XDIAG_THROW("Invalid stride given. This is a bug, please report.");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// This initialization copies the memory (complex)
void State::initcopy(const complex *ptr, int64_t nrows, int64_t ncols) try {
  init0(false, nrows, ncols);
  std::copy(ptr, ptr + size(), reinterpret_cast<complex *>(storage_.data()));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, bool real, int64_t ncols) try
    : valid_(true), block_(block) {
  init0(real && xdiag::isreal(block), xdiag::size(block), ncols);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, arma::vec const &vector) try
    : valid_(true), block_(block) {
  if (vector.size() != xdiag::size(block)) {
    XDIAG_THROW(
        "Size of block does not agree with size of given armadillo vector");
  }
  if (xdiag::isreal(block)) {
    initcopy(vector.memptr(), vector.size(), 1); // real init
  } else {
    arma::cx_vec vectorc = utils::to_cx_vec(vector);
    initcopy(vectorc.memptr(), vector.size(), 1); // cplx init
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, arma::cx_vec const &vector) try
    : valid_(true), block_(block) {
  if (vector.size() != xdiag::size(block)) {
    XDIAG_THROW(
        "Size of block does not agree with size of given armadillo vector");
  }
  initcopy(vector.memptr(), vector.size(), 1); // cplx init
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, arma::mat const &matrix) try
    : valid_(true), block_(block) {
  if (matrix.n_rows != xdiag::size(block)) {
    XDIAG_THROW("Size of block does not agree with number of rows of given "
                "armadillo matrix");
  }
  if (xdiag::isreal(block)) {
    initcopy(matrix.memptr(), matrix.n_rows, matrix.n_cols); // real init
  } else {
    arma::cx_mat matrixc = utils::to_cx_mat(matrix);
    initcopy(matrixc.memptr(), matrix.n_rows, matrix.n_cols); // cplx init
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, arma::cx_mat const &matrix) try
    : valid_(true), block_(block) {
  if (matrix.n_rows != xdiag::size(block)) {
    XDIAG_THROW("Size of block does not agree with number of rows of given "
                "armadillo matrix");
  }
  initcopy(matrix.memptr(), matrix.n_rows, matrix.n_cols); // cplx init
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, double const *ptr, int64_t ncols,
             int64_t stride) try
    : valid_(true), block_(block) {
  initcopy(ptr, xdiag::size(block), ncols, stride); // real init
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

State::State(Block const &block, complex const *ptr, int64_t ncols) try
    : valid_(true), block_(block) {
  initcopy(ptr, xdiag::size(block), ncols); // complex init
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool State::isvalid() const { return valid_; }
int64_t State::nsites() const { return xdiag::nsites(block_); }
bool State::isreal() const { return real_; }

State State::real() const try {
  if (isreal()) {
    return (*this);
  } else {
    double *ptr = storage_.data();
    return std::visit(
        [&](auto &&block) { return State(block, ptr, ncols_, 2); }, block_);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return (*this);
}

State State::imag() const try { // TODO: DOES THIS DO WHAT IT"S SUPPOSED TO?
  if (isreal()) {
    return State(block_, true, ncols_);
  } else {
    double *ptr = storage_.data();
    return std::visit(
        [&](auto &&block) { return State(block, ptr + 1, ncols_, 2); }, block_);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return (*this);
}

void State::make_complex() try {
  if (isreal()) {
    real_ = false;

    safe_resize(storage_, 2 * size());

    double *ptr = storage_.data();
    for (int64_t i = size() - 1; i >= 0; --i) {
      std::swap(ptr[i << 1], ptr[i]);
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t State::dim() const { return xdiag::dim(block_); }

int64_t State::size() const { return nrows_ * ncols_; }
int64_t State::nrows() const { return nrows_; }
int64_t State::ncols() const { return ncols_; }
Block State::block() const { return block_; }

State State::col(int64_t n, bool copy) const try {
  if (n >= ncols_) {
    XDIAG_THROW("Column index larger than the number of columns");
  }
  if (real_) {
    return std::visit([&](auto &&b) { return State(b, vector(n, copy)); },
                      block_);
  } else {
    return std::visit([&](auto &&b) { return State(b, vectorC(n, copy)); },
                      block_);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return State();
}

arma::vec State::vector(int64_t n, bool copy) const try {
  if (!real_) {
    XDIAG_THROW("Cannot return a real armadillo vector from a "
                "complex State (maybe use vectorC(...) instead)");
  } else if (n >= ncols_) {
    XDIAG_THROW("Column index larger than the number of columns");
  } else if (n < 0) {
    XDIAG_THROW("Negative column index");
  }
  return arma::vec(storage_.data() + n * nrows_, nrows_, copy, !copy);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return arma::vec();
}

arma::mat State::matrix(bool copy) const try {
  if (!real_) {
    XDIAG_THROW("Cannot return a real armadillo matrix from a "
                "complex state (maybe use matrixC(...) instead)");
  }
  return arma::mat(storage_.data(), nrows_, ncols_, copy, !copy);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return arma::mat();
}

arma::cx_vec State::vectorC(int64_t n, bool copy) const try {
  if (real_) {
    XDIAG_THROW("Cannot return a complex armadillo vector from a "
                "real State (maybe use vector(...) instead)");
  } else if (n >= ncols_) {
    XDIAG_THROW("Column index larger than the number of columns");
  } else if (n < 0) {
    XDIAG_THROW("Negative column index");
  }

  return arma::cx_vec(reinterpret_cast<complex *>(storage_.data()) + n * nrows_,
                      nrows_, copy, !copy);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return arma::cx_vec();
}

arma::cx_mat State::matrixC(bool copy) const try {
  if (real_) {
    XDIAG_THROW("Cannot return a complex armadillo matrix from a "
                "real state (maybe use matrix(...) instead)");
  }
  return arma::cx_mat(reinterpret_cast<complex *>(storage_.data()), nrows_,
                      ncols_, copy, !copy);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
  return arma::cx_mat();
}

double *State::memptr() { return storage_.data(); }
complex *State::memptrC() {
  return reinterpret_cast<complex *>(storage_.data());
}
double *State::colptr(int64_t col) {
  if ((col < 0) || (col >= ncols_)) {
    XDIAG_THROW("Invalid column index requested");
  }
  return memptr() + col * nrows_;
}
complex *State::colptrC(int64_t col) {
  if ((col < 0) || (col >= ncols_)) {
    XDIAG_THROW("Invalid column index requested");
  }
  return memptrC() + col * nrows_;
}

bool isvalid(State const &s) { return s.isvalid(); }
int64_t nsites(State const &s) { return s.nsites(); }
bool isapprox(State const &v, State const &w, double rtol, double atol) try {
  if (v.block() == w.block()) {
    if (isreal(v) && isreal(w)) {
      return arma::approx_equal(v.matrix(false), w.matrix(false), "both", atol,
                                rtol);
    } else if (!isreal(v) && !isreal(w)) {
      return arma::approx_equal(v.matrixC(false), w.matrixC(false), "both",
                                atol, rtol);
    } else {
      return false;
    }
  } else {
    return false;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

bool isreal(State const &s) { return s.isreal(); }
State real(State const &s) { return s.real(); }
State imag(State const &s) { return s.imag(); }
void make_complex(State &s) { return s.make_complex(); }
int64_t dim(State const &s) { return s.dim(); }
int64_t size(State const &s) { return s.size(); }
int64_t nrows(State const &s) { return s.nrows(); }
int64_t ncols(State const &s) { return s.ncols(); }

State col(State const &s, int64_t n, bool copy) { return s.col(n, copy); }
arma::vec vector(State const &s, int64_t n, bool copy) {
  return s.vector(n, copy);
}
arma::mat matrix(State const &s, bool copy) { return s.matrix(copy); }
arma::cx_vec vectorC(State const &s, int64_t n, bool copy) {
  return s.vectorC(n, copy);
}
arma::cx_mat matrixC(State const &s, bool copy) { return s.matrixC(copy); }

std::ostream &operator<<(std::ostream &out, State const &state) {
  if (state.isvalid()) {
    if (state.isreal()) {
      out << "REAL State\n";
    } else {
      out << "COMPLEX State\n";
    }
    out << "Block:\n";
    out << state.block();
  } else {
    out << "INVALID State\n";
  }
  return out;
}
std::string to_string(State const &state) { return to_string_generic(state); }

State &operator*=(State &X, double alpha) try {
  if (isreal(X)) {
    X.matrix(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(State const &X, double alpha) try {
  auto res = X;
  res *= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(double alpha, State const &X) try {
  return X * alpha;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator*=(State &X, complex alpha) try {
  if (isreal(X)) {
    X.make_complex();
    X.matrixC(false) *= alpha;
  } else {
    X.matrixC(false) *= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(State const &X, complex alpha) try {
  auto res = X;
  res *= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator*(complex alpha, State const &X) try {
  return X * alpha;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator/=(State &X, double alpha) try {
  if (isreal(X)) {
    X.matrix(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}
State operator/(State const &X, double alpha) try {
  auto res = X;
  res /= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator/=(State &X, complex alpha) try {
  if (isreal(X)) {
    X.make_complex();
    X.matrixC(false) /= alpha;
  } else {
    X.matrixC(false) /= alpha;
  }
  return X;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator/(State const &X, complex alpha) try {
  auto res = X;
  res /= alpha;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator+=(State &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot add two states from different blocks");
  }
  if (v.ncols() != w.ncols()) {
    XDIAG_THROW("Cannot add two states with different number of columns");
  }

  if (isreal(v) && isreal(w)) {
    v.matrix(false) += w.matrix(false);
  } else if (!isreal(v) && isreal(w)) {
    auto wcplx = w;
    wcplx.make_complex();
    v.matrixC(false) += wcplx.matrixC(false);
  } else if (isreal(v) && !isreal(w)) {
    v.make_complex();
    v.matrixC(false) += w.matrixC(false);
  } else { // (!isreal(v) && !isreal(w))
    v.matrixC(false) += w.matrixC(false);
  }

  return v;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State &operator-=(State &v, State const &w) try {
  if (v.block() != w.block()) {
    XDIAG_THROW("Cannot subtract two states from different blocks");
  }
  if (v.ncols() != w.ncols()) {
    XDIAG_THROW("Cannot subtract two states with different number of columns");
  }

  if (isreal(v) && isreal(w)) {
    v.matrix(false) -= w.matrix(false);
  } else if (!isreal(v) && isreal(w)) {
    auto wcplx = w;
    wcplx.make_complex();
    v.matrixC(false) -= wcplx.matrixC(false);
  } else if (isreal(v) && !isreal(w)) {
    v.make_complex();
    v.matrixC(false) -= w.matrixC(false);
  } else { // (!isreal(v) && !isreal(w))
    v.matrixC(false) -= w.matrixC(false);
  }

  return v;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator+(State const &v, State const &w) try {
  auto res = v;
  res += w;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator-(State const &v, State const &w) try {
  auto res = v;
  res -= w;
  return res;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

State operator-(State const &v) try {
  return (-1.0) * v;
} catch (Error const &error) {
  XDIAG_RETHROW(error);
}

} // namespace xdiag
