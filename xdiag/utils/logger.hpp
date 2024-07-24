#pragma once

#ifdef XDIAG_USE_MPI
#include <xdiag/utils/logger_mpi.hpp>
#else
#include <xdiag/utils/logger_serial.hpp>
#endif

namespace xdiag {

#ifdef XDIAG_USE_MPI
inline auto & Log = LogMPI;
#else
inline auto & Log = LogSerial;
#endif

void set_verbosity(int64_t level);
  
}
