#pragma once
#ifdef XDIAG_USE_MPI

#include <mpi.h>
#include <xdiag/parallel/mpi/logger_mpi.h>

namespace xdiag {

double rightnow_mpi();
void timing_mpi(double const &t0, double const &t1, std::string msg = "",
                int verbosity = 1);
void tic_mpi(bool begin = true, std::string msg = "", int verbosity = 1);
void toc_mpi(std::string msg = "", int verbosity = 1);

} // namespace xdiag
#endif
