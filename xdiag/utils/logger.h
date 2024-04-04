#pragma once

#ifdef XDIAG_USE_MPI
#include <xdiag/parallel/mpi/logger_mpi.h>
#else
#include <xdiag/utils/logger_serial.h>
#endif

namespace xdiag {

#ifdef XDIAG_USE_MPI
inline auto & Log = LogMPI;
#else
inline auto & Log = LogSerial;
#endif

}
