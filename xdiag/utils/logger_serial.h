#pragma once

#include <cstdlib>
#include <iostream>

#define FMT_HEADER_ONLY
#include <xdiag/extern/fmt/format.h>

#include <xdiag/utils/error.h>

namespace xdiag {

class Logger {
public:
  Logger() : verbosity_(0){};

  inline void set_verbosity(int verbosity) { verbosity_ = verbosity; }
  inline int verbosity() { return verbosity_; }

  template <typename... Args>
  inline void out(const std::string &format, const Args &...args) try {
    std::cout << fmt::format(format, args...) << "\n";
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void warn(const std::string &format, const Args &...args) try {
    std::cout << fmt::format(format, args...) << "\n" << std::flush;
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void err(const std::string &format, const Args &...args) try {
    std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void out(int level, const std::string &format, const Args &...args) try {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n";
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void warn(int level, const std::string &format, const Args &...args) try {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void err(int level, const std::string &format, const Args &...args) try {
    if (level <= verbosity_)
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void operator()(const std::string &format, const Args &...args) try {
    out(format, args...);
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

  template <typename... Args>
  inline void operator()(int level, const std::string &format,
                  const Args &...args) try {
    out(level, format, args...);
  } catch (...) {
    XDiagRethrow("Unable to print output using Logger");
  }

private:
  int verbosity_;
};

inline Logger LogSerial;

} // namespace xdiag