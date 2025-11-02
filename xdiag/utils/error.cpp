// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "error.hpp"
#ifdef XDIAG_USE_MPI
#include <mpi.h>
#endif
#include <iomanip> // std::setw
#include <variant>

#define FMT_HEADER_ONLY
#include <xdiag/extern/fmt/format.hpp>

#ifndef XDIAG_DISABLE_COLOR
#include <xdiag/extern/fmt/color.hpp>
#endif

#include <xdiag/extern/armadillo/armadillo>

#ifdef __APPLE__
#ifdef __clang__
/// DIRTY HACK to make std::visit work on old MacOS versions
const char *std::bad_variant_access::what() const noexcept {
  return "bad_variant_access";
}
#endif
#endif

namespace xdiag {

Error::Error(std::string message, std::string original_message)
    : original_message_(original_message), messages_({message}) {
  if (original_message == "") {
    original_message_ = message;
  }
}
Error::Error(Error const &error, std::string message)
    : original_message_(error.original_message()), messages_(error.messages()) {
  messages_.push_back(message);
}
// Error::~Error(){};

std::vector<std::string> const &Error::messages() const { return messages_; }
std::string Error::original_message() const { return original_message_; }

const char *Error::what() const noexcept { return original_message_.c_str(); }

std::string cut_file(const char *file) {
  std::string dir(XDIAG_DIRECTORY);
  int n = dir.size();
  return std::string(file + n + 1);
}

void throw_error(std::string message, const char *file, const char *func,
                 int line) {
  std::ostringstream ss;
  ss << fmt::format(fmt::emphasis::bold, func) << " (" << cut_file(file) << ":"
     << line << ")";
  throw Error(ss.str(), message);
}

void rethrow_error(Error const &error, const char *file, const char *func,
                   int line) {

#ifdef XDIAG_DISABLE_COLOR
  std::string funcname(func);
#else
  std::string funcname = fmt::format(fmt::emphasis::bold, func);
#endif

  std::ostringstream ss;
  ss << fmt::format(fmt::emphasis::bold, func) << " (" << cut_file(file) << ":"
     << line << ")";
  throw Error(error, ss.str());
}

void error_trace(Error const &error) {

#ifdef XDIAG_USE_MPI
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank == 0) {
#endif

#ifdef XDIAG_DISABLE_COLOR
    std::string error_head = "XDiag ERROR: ";
#else
  std::string error_head = fmt::format(
      fg(fmt::color::crimson) | fmt::emphasis::bold, "XDiag ERROR: ");
#endif

    std::cerr << error_head << error.what() << "\n";
    std::cerr << "Stacktrace (C++):\n";
    int64_t idx = 1;
    for (auto message : error.messages()) {
      if (idx < 10) {
        std::cerr << " [" << idx << "] " << message << "\n";
      } else {
        std::cerr << "[" << idx << "] " << message << "\n";
      }
      ++idx;
    }

#ifdef XDIAG_USE_MPI
  }
#endif
}

void check_dimension_works_with_blas_int_size(int64_t dim) try {
  // Backend 32 bit Blas implementation
  if ((sizeof(arma::blas_int) == 4) && (dim > ((int64_t)1 << 31) - 2)) {
    XDIAG_THROW("Trying to create a block whose dimension is too large for the "
                "backend BLAS routines. Your backend BLAS routines have a 32 "
                "bit interface, and a vector beyond a dimension of 2^31 cannot "
                "be represented. If you want to still perform this "
                "calculation, you will need to compile XDiag with a 64 bit "
                "backend, for example the ILP64 interface of IntelMKL.");
  }
  // Backend 64 bit Blas implementation
  else if ((sizeof(arma::blas_int) == 8) && (dim > ((int64_t)1 << 62) - 2)) {
    XDIAG_THROW(
        "Trying to create a block whose dimension is too large for the "
        "backend BLAS routines. The block dimension requested is larger than"
        "2^63 which anyway is not feasible");
  }
}
XDIAG_CATCH

template <typename bit_t> void check_nsites_work_with_bits(int64_t nsites) try {
  int64_t n_bits = std::numeric_limits<bit_t>::digits;
  if (nsites >= n_bits) {
    XDIAG_THROW(
        fmt::format("Cannot encode basis with nsites={} using only {} bits. "
                    "Consider using a different backend if possible.",
                    nsites, n_bits));
  }
}
XDIAG_CATCH

template void check_nsites_work_with_bits<uint16_t>(int64_t nsites);
template void check_nsites_work_with_bits<uint32_t>(int64_t nsites);
template void check_nsites_work_with_bits<uint64_t>(int64_t nsites);
} // namespace xdiag
