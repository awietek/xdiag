#pragma once

#include <iostream>

#ifdef XDIAG_USE_MPI
#include <mpi.h>
#endif

namespace xdiag::utils {

template <typename T> inline void print_pretty(const char *id, T X) {
  std::cout << id << ":\n";
  std::cout << X << "\n";
}

#ifdef XDIAG_USE_MPI
template <typename T> inline void print_pretty_mpi(const char *id, T X) {
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    print_pretty(id, X);
  }
}
#endif

} // namespace xdiag::utils

// XDIAG_SHOW belongs to XDIAG_API
#ifdef XDIAG_USE_MPI
#define XDIAG_SHOW(X) xdiag::utils::print_pretty_mpi(#X, X)
#else
#define XDIAG_SHOW(X) xdiag::utils::print_pretty(#X, X)
#endif
