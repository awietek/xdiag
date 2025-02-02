#pragma once
#include <chrono>
#include <string>

#include <xdiag/utils/logger.hpp>

namespace xdiag {

using namespace std::chrono;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;

auto inline rightnow() -> decltype(high_resolution_clock::now()) {
  return high_resolution_clock::now();
}

template <class Clock, class Duration = typename Clock::duration>
inline void timing(time_point<Clock, Duration> const &t0,
                   time_point<Clock, Duration> const &t1, std::string msg = "",
                   int verbosity = 0) {
  auto td = duration_cast<microseconds>(t1 - t0).count();
  double tds = (double)td / 1000000;
  if (msg != "") {
    Log.out(verbosity, "{}: {:.5f} secs", msg, tds);
  } else {
    Log.out(verbosity, "{:.5f} secs", tds);
  }
}

inline void tic(bool begin = true, std::string msg = "", int verbosity = 0) {
  static auto t0 = rightnow();
  if (begin) {
    t0 = rightnow();
  } else {
    auto t1 = rightnow();
    timing(t0, t1, msg, verbosity);
  }
}

inline void toc(std::string msg = "", int verbosity = 0) {
  tic(false, msg, verbosity);
}

} // namespace xdiag
