#include "fill.hpp"

#include <variant>

#include <xdiag/blocks/electron/electron.hpp>
#include <xdiag/blocks/spinhalf/spinhalf.hpp>
#include <xdiag/blocks/tj/tj.hpp>

namespace xdiag {

template <typename basis_t, typename coeff_t, class fill_f>
void fill_generic(basis_t const &basis, coeff_t *ptr, fill_f func) try {
  int64_t idx = 0;
  for (auto bits : basis) {
    ptr[idx++] = func(bits);
  }
  assert(idx == basis.size());
} catch (...) {
  rethrow(__func__,
          std::runtime_error("Unable to perform generic fill operation"));
}

void fill(State &state, std::function<double(uint64_t)> func, int64_t col) try {
  auto block = state.block();

  // Check whether one parameter function is warranted
  if (std::holds_alternative<tJ>(block)) {
    throw(std::logic_error(
        "Cannot fill tJ block with function of one argument only. The function "
        "needs to accept an std::pair<bit_t, bit_t> instead."));
  } else if (std::holds_alternative<Electron>(block)) {
    throw(std::logic_error(
        "Cannot fill Electron block with function of one "
        "argument only. The function "
        "needs to accept an std::pair<bit_t, bit_t> instead."));
  }

  auto basis = block.basis();

  if (state.isreal()) {
    double *colptr = state.colptr(col);
    fill_generic(basis, colptr, func);
  } else {
    complex *colptr = state.colptrC(col);
    fill_generic(basis, colptr, func);
  }
} catch (...) {
  rethrow(__func__, std::runtime_error("Unable to fill State"));
}

void fill(State &state,
          std::function<double(std::pair<uint64_t, uint64_t>)> func,
          int64_t col) try {
  auto block = state.block();

  // Check whether one parameter function is warranted
  if (std::holds_alternative<Spinhalf>(Spinhalf)) {
    throw(std::logic_error("Cannot fill Spinhalf block with function of "
                           "std::pair<bit_t, bit_t>."));
  }

  auto basis = block.basis();

  if (state.isreal()) {
    double *colptr = state.colptr(col);
    fill_generic(basis, colptr, func);
  } else {
    complex *colptr = state.colptrC(col);
    fill_generic(basis, colptr, func);
  }
} catch (...) {
  rethrow(__func__, std::runtime_error("Unable to fill State"));
}

void fill(State &state, std::function<complex(uint64_t)> func,
          int64_t col) try {
  auto block = state.block();

  // Check whether one parameter function is warranted
  if (std::holds_alternative<tJ>(block)) {
    throw(std::logic_error(
        "Cannot fill tJ block with function of one argument only. The function "
        "needs to accept an std::pair<bit_t, bit_t> instead."));
  } else if (std::holds_alternative<Electron>(block)) {
    throw(std::logic_error(
        "Cannot fill Electron block with function of one "
        "argument only. The function "
        "needs to accept an std::pair<bit_t, bit_t> instead."));
  }

  auto basis = block.basis();

  if (state.isreal()) {
    throw("Cannot fill real state with a function that returns complex values");
  } else {
    complex *colptr = state.colptrC(col);
    fill_generic(basis, colptr, func);
  }
} catch (...) {
  rethrow(__func__, std::runtime_error("Unable to fill State"));
}

void fill(State &state,
          std::function<complex(std::pair<uint64_t, uint64_t>)> func,
          int64_t col) try {
  auto block = state.block();

  // Check whether one parameter function is warranted
  if (std::holds_alternative<Spinhalf>(Spinhalf)) {
    throw(std::logic_error("Cannot fill Spinhalf block with function of "
                           "std::pair<bit_t, bit_t>."));
  }

  auto basis = block.basis();

  if (state.isreal()) {
    throw("Cannot fill real state with a function that returns complex values");
  } else {
    complex *colptr = state.colptrC(col);
    fill_generic(basis, colptr, func);
  }
} catch (...) {
  rethrow(__func__, std::runtime_error("Unable to fill State"));
}

} // namespace xdiag
