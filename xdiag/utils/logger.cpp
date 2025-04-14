// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "logger.hpp"

namespace xdiag {

void set_verbosity(int64_t level) { xdiag::Log.set_verbosity(level); }

} // namespace xdiag
