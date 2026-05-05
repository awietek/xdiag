// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "timing.hpp"

#include <xdiag/utils/logger.hpp>

namespace xdiag {

auto rightnow() -> decltype(clock_t::now()) { return clock_t::now(); }

void timing(std::chrono::time_point<clock_t> const &t0,
            std::chrono::time_point<clock_t> const &t1, std::string msg,
            int verbosity) {
  auto td =
      std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
  double tds = (double)td / 1000000;
  if (msg != "") {
    Log(verbosity, "{}: {:.5f} secs", msg, tds);
  } else {
    Log(verbosity, "{:.5f} secs", tds);
  }
}

void tic(bool begin, std::string msg, int verbosity) {
  static auto t0 = rightnow();
  if (begin) {
    t0 = rightnow();
  } else {
    auto t1 = rightnow();
    timing(t0, t1, msg, verbosity);
  }
}

void toc(std::string msg, int verbosity) { tic(false, msg, verbosity); }

} // namespace xdiag
