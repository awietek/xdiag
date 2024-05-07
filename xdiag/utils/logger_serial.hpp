#pragma once

#include <cstdlib>
#include <iostream>

#define FMT_HEADER_ONLY
#include <xdiag/extern/fmt/format.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag {

class Logger {
public:
  Logger() : verbosity_(0){};

  inline void set_verbosity(int verbosity) { verbosity_ = verbosity; }
  inline int verbosity() { return verbosity_; }

  template <typename... Args>
  inline void out(const std::string &format, const Args &...args) {
    std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  inline void warn(const std::string &format, const Args &...args) {
    std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  inline void err(const std::string &format, const Args &...args) {
    std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  }

  template <typename... Args>
  inline void out(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  inline void warn(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  inline void err(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  }

  template <typename... Args>
  inline void operator()(const std::string &format, const Args &...args) {
    out(format, args...);
  }

  template <typename... Args>
  inline void operator()(int level, const std::string &format,
                         const Args &...args) {
    out(level, format, args...);
  }

private:
  int verbosity_;
};

inline Logger LogSerial;

} // namespace xdiag
