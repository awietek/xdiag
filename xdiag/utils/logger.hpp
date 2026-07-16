// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#ifdef XDIAG_DISTRIBUTED
#include <xdiag/utils/logger_mpi.hpp>
#else
#include <xdiag/utils/logger_serial.hpp>
#endif

#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

#ifdef XDIAG_DISTRIBUTED
XDIAG_API inline auto &Log = LogMPI;
#else
XDIAG_API inline auto &Log = LogSerial;
#endif

XDIAG_API void set_verbosity(int64_t level);

} // namespace xdiag
