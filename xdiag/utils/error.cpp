#include "error.hpp"

#include <iomanip> // std::setw
#include <variant>

#ifndef XDIAG_DISABLE_COLOR
#include <xdiag/extern/fmt/color.hpp>
#endif

#ifdef __APPLE__
#ifdef __clang__
/// DIRTY HACK to make std::visit work on old MacOS versions
const char *std::bad_variant_access::what() const noexcept {
  return "bad_variant_access";
}
#endif
#endif

namespace xdiag {

Error::Error(std::string message, std::string original_message)
    : original_message_(original_message), messages_({message}) {
  if (original_message == "") {
    original_message_ = message;
  }
}
Error::Error(Error const &error, std::string message)
    : original_message_(error.original_message()), messages_(error.messages()) {
  messages_.push_back(message);
}

std::vector<std::string> const &Error::messages() const { return messages_; }
std::string Error::original_message() const { return original_message_; }

const char *Error::what() const noexcept { return original_message_.c_str(); }

std::string cut_file(const char *file) {
  std::string dir(XDIAG_DIRECTORY);
  int n = dir.size();
  return std::string(file + n + 1);
}

void throw_error(std::string message, const char *file, const char *func,
                 int line) {
  std::ostringstream ss;
  ss << fmt::format(fmt::emphasis::bold, func) << " (" << cut_file(file) << ":"
     << line << ")";
  throw Error(ss.str(), message);
}

void rethrow_error(Error const &error, const char *file, const char *func,
                   int line) {

#ifdef XDIAG_DISABLE_COLOR
  std::string funcname(func);
#else
  std::string funcname = fmt::format(fmt::emphasis::bold, func);
#endif

  std::ostringstream ss;
  ss << fmt::format(fmt::emphasis::bold, func) << " (" << cut_file(file) << ":"
     << line << ")";
  throw Error(error, ss.str());
}

void error_trace(Error const &error) {

#ifdef XDIAG_DISABLE_COLOR
  std::string error_hear = "XDiag ERROR: ";
#else
  std::string error_head = fmt::format(
      fg(fmt::color::crimson) | fmt::emphasis::bold, "XDiag ERROR: ");
#endif

  std::cerr << error_head << error.what() << "\n";
  std::cerr << "Stacktrace (C++):\n";
  int64_t idx = 1;
  for (auto message : error.messages()) {
    if (idx < 10) {
      std::cerr << " [" << idx << "] " << message << "\n";
    } else {
      std::cerr << "[" << idx << "] " << message << "\n";

    }
    ++idx;
  }
}

} // namespace xdiag
