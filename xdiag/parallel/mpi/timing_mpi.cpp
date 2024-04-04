#include "timing_mpi.hpp"

#include <mpi.h>

namespace xdiag {

double rightnow_mpi() { return MPI_Wtime(); }

void timing_mpi(double const &t0, double const &t1, std::string msg,
                int verbosity) {
  auto td = t1 - t0;
  if (msg != "")
    LogMPI.out(verbosity, "{}: {:.6f} secs", msg, td);
  else
    LogMPI.out(verbosity, "{:.6f} secs", td);
}

void tic_mpi(bool begin, std::string msg, int verbosity) {
  static double t0;
  if (begin) {
    t0 = rightnow_mpi();
  } else {
    auto t1 = rightnow_mpi();
    timing_mpi(t0, t1, msg, verbosity);
  }
}

void toc_mpi(std::string msg, int verbosity) { tic_mpi(false, msg, verbosity); }

} // namespace xdiag
