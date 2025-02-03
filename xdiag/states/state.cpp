#include "state.hpp"

#include <xdiag/blocks/blocks.hpp>

namespace xdiag {

static void safe_resize(std::vector<double> &vec, int64_t size) try {
  vec.resize(size);
} catch (...) {
  XDIAG_THROW("Unable to resize vector");
}

State::State(Block const &block, bool real, int64_t ncols) try
    : real_(real), dim_(xdiag::dim(block)), nrows_(xdiag::size(block)),
      ncols_(ncols), block_(block) {
  if (real) {
    safe_resize(storage_, size());
  } else {
    safe_resize(storage_, 2 * size());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
State::State(block_t const &block, bool real, int64_t ncols) try
    : real_(real), dim_(block.dim()), nrows_(block.size()), ncols_(ncols),
      block_(block) {
  if (real) {
    safe_resize(storage_, size());
  } else {
    safe_resize(storage_, 2 * size());
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
State::State(block_t const &block, double const *ptr, int64_t ncols,
             int64_t stride) try
    : real_(true), nrows_(block.size()), ncols_(ncols), block_(block) {

  safe_resize(storage_, size());

  if (stride == 1) {
    std::copy(ptr, ptr + size(), storage_.data());
  } else {
    for (int64_t i = 0, is = 0; i < size(); ++i, is += stride) {
      storage_[i] = ptr[is];
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
State::State(block_t const &block, complex const *ptr, int64_t ncols) try
    : real_(false), nrows_(block.size()), ncols_(ncols), block_(block) {

  safe_resize(storage_, 2 * size());
  std::copy(ptr, ptr + size(), reinterpret_cast<complex *>(storage_.data()));

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t, typename coeff_t>
State::State(block_t const &block, arma::Col<coeff_t> const &vector) try
    : real_(xdiag::isreal<coeff_t>()), nrows_(block.size()), ncols_(1),
      block_(block) {
  if (block.size() != (int64_t)vector.n_rows) {
    XDIAG_THROW("Block dimension not equal to vector dimension");
  }

  if (real_) {
    safe_resize(storage_, size());
  } else {
    safe_resize(storage_, 2 * size());
  }

  std::copy(vector.memptr(), vector.memptr() + size(),
            reinterpret_cast<coeff_t *>(storage_.data()));

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t, typename coeff_t>
State::State(block_t const &block, arma::Mat<coeff_t> const &matrix) try
    : real_(xdiag::isreal<coeff_t>()), nrows_(matrix.n_rows),
      ncols_(matrix.n_cols), block_(block) {

  if (block.size() != (int64_t)matrix.n_rows) {
    XDIAG_THROW("Block dimension not equal to number of rows in matrix");
  }

  if (real_) {
    safe_resize(storage_, size());
  } else {
    safe_resize(storage_, 2 * size());
  }

  std::copy(matrix.memptr(), matrix.memptr() + size(),
            reinterpret_cast<coeff_t *>(storage_.data()));

} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

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

template XDIAG_API State::State(Spinhalf const &, bool, int64_t);
template XDIAG_API State::State(tJ const &, bool, int64_t);
template XDIAG_API State::State(Electron const &, bool, int64_t);

template XDIAG_API State::State(Spinhalf const &, double const *, int64_t,
                                int64_t);
template XDIAG_API State::State(tJ const &, double const *, int64_t, int64_t);
template XDIAG_API State::State(Electron const &, double const *, int64_t,
                                int64_t);

template XDIAG_API State::State(Spinhalf const &, complex const *, int64_t);
template XDIAG_API State::State(tJ const &, complex const *, int64_t);
template XDIAG_API State::State(Electron const &, complex const *, int64_t);

template XDIAG_API State::State(Spinhalf const &, arma::Col<double> const &);
template XDIAG_API State::State(tJ const &, arma::Col<double> const &);
template XDIAG_API State::State(Electron const &, arma::Col<double> const &);
template XDIAG_API State::State(Spinhalf const &, arma::Col<complex> const &);
template XDIAG_API State::State(tJ const &, arma::Col<complex> const &);
template XDIAG_API State::State(Electron const &, arma::Col<complex> const &);

template XDIAG_API State::State(Spinhalf const &, arma::Mat<double> const &);
template XDIAG_API State::State(tJ const &, arma::Mat<double> const &);
template XDIAG_API State::State(Electron const &, arma::Mat<double> const &);
template XDIAG_API State::State(Spinhalf const &, arma::Mat<complex> const &);
template XDIAG_API State::State(tJ const &, arma::Mat<complex> const &);
template XDIAG_API State::State(Electron const &, arma::Mat<complex> const &);

#ifdef XDIAG_USE_MPI
template XDIAG_API State::State(SpinhalfDistributed const &, bool, int64_t);
template XDIAG_API State::State(SpinhalfDistributed const &, double const *,
                                int64_t, int64_t);
template XDIAG_API State::State(SpinhalfDistributed const &, complex const *,
                                int64_t);
template XDIAG_API State::State(SpinhalfDistributed const &,
                                arma::Col<double> const &);
template XDIAG_API State::State(SpinhalfDistributed const &,
                                arma::Col<complex> const &);
template XDIAG_API State::State(SpinhalfDistributed const &,
                                arma::Mat<double> const &);
template XDIAG_API State::State(SpinhalfDistributed const &,
                                arma::Mat<complex> const &);

template XDIAG_API State::State(tJDistributed const &, bool, int64_t);
template XDIAG_API State::State(tJDistributed const &, double const *, int64_t,
                                int64_t);
template XDIAG_API State::State(tJDistributed const &, complex const *,
                                int64_t);
template XDIAG_API State::State(tJDistributed const &,
                                arma::Col<double> const &);
template XDIAG_API State::State(tJDistributed const &,
                                arma::Col<complex> const &);
template XDIAG_API State::State(tJDistributed const &,
                                arma::Mat<double> const &);
template XDIAG_API State::State(tJDistributed const &,
                                arma::Mat<complex> const &);

#endif

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
  if (state.isreal()) {
    out << "REAL State\n";
  } else {
    out << "COMPLEX State\n";
  }
  out << "Block:\n";
  out << state.block();
  return out;
}
std::string to_string(State const &state) { return to_string_generic(state); }
} // namespace xdiag
