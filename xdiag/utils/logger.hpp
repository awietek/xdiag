#pragma once

#include <xdiag/common.hpp>

#ifdef XDIAG_USE_MPI
#include <xdiag/utils/logger_mpi.hpp>
#else
#include <xdiag/utils/logger_serial.hpp>
#endif

namespace xdiag {

#ifdef XDIAG_USE_MPI
XDIAG_API inline auto &Log = LogMPI;
#else
XDIAG_API inline auto &Log = LogSerial;
#endif

XDIAG_API void set_verbosity(int64_t level);

} // namespace xdiag

