#pragma once

#include <cstdlib>
#include <iostream>

#define FMT_HEADER_ONLY
#include <extern/fmt/format.h>

namespace hydra {

class Logger {
public:
  Logger() : verbosity_(0){};

  void set_verbosity(int verbosity) { verbosity_ = verbosity; }
  int verbosity() { return verbosity_; }

  template <typename... Args>
  void out(const std::string &format, const Args &...args) {
    std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  void warn(const std::string &format, const Args &...args) {
    std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  void err(const std::string &format, const Args &...args) {
    std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  }

  template <typename... Args>
  void out(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n";
  }

  template <typename... Args>
  void warn(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cout << fmt::format(format, args...) << "\n" << std::flush;
  }

  template <typename... Args>
  void err(int level, const std::string &format, const Args &...args) {
    if (level <= verbosity_)
      std::cerr << fmt::format(format, args...) << "\n" << std::flush;
    exit(EXIT_FAILURE);
  }

  template <typename... Args>
  void operator()(const std::string &format, const Args &...args) {
    out(format, args...);
  }

  template <typename... Args>
  void operator()(int level, const std::string &format, const Args &...args) {
    out(level, format, args...);
  }

private:
  int verbosity_;
};

inline Logger LogSerial;

} // namespace hydra
