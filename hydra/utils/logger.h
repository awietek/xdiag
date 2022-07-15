#pragma once

#ifdef HYDRA_ENABLE_MPI
#include <mpi.h>
#include <hydra/mpi/logger_mpi.h>
#else
#include <lila/utils/logger.h>
#endif

namespace hydra {

#ifdef HYDRA_ENABLE_MPI
inline auto & Log = LogMPI;
#else
inline auto & Log = lila::LogSerial;
#endif

}
