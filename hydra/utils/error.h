#pragma once

#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace hydra {

template <class... Args>
inline void rethrow2(const char *file, const char *func, int line,
                     Args &&...args) {
  // build an error message
  std::ostringstream ss;
  ss << "File    : " << file << "\n";
  ss << "       Function: " << func << "\n";
  ss << "       Line    : " << line << "\n";
  ss << "       ";
  auto sep = " : ";
  using expand = int[];
  void(expand{0, ((ss << args), sep = ", ", 0)...});
  // figure out what kind of exception is active
  try {
    std::rethrow_exception(std::current_exception());
  } catch (const std::invalid_argument &e) {
    std::throw_with_nested(std::invalid_argument(ss.str()));
  } catch (const std::logic_error &e) {
    std::throw_with_nested(std::logic_error(ss.str()));
  }
  // etc - default to a runtime_error
  catch (...) {
    std::throw_with_nested(std::runtime_error(ss.str()));
  }
}

template <class Context, class... Args>
inline void rethrow(Context &&context, Args &&...args) {
  // build an error message
  std::ostringstream ss;
  ss << context;
  auto sep = " : ";
  using expand = int[];
  void(expand{0, ((ss << sep << args), sep = ", ", 0)...});
  // figure out what kind of exception is active
  try {
    std::rethrow_exception(std::current_exception());
  } catch (const std::invalid_argument &e) {
    std::throw_with_nested(std::invalid_argument(ss.str()));
  } catch (const std::logic_error &e) {
    std::throw_with_nested(std::logic_error(ss.str()));
  }
  // etc - default to a runtime_error
  catch (...) {
    std::throw_with_nested(std::runtime_error(ss.str()));
  }
}

// unwrap nested exceptions, printing each nested exception to
// std::cerr
void traceback(const std::exception &e, std::size_t depth = 0);

} // namespace hydra

#define HydraRethrow(msg) hydra::rethrow2(__FILE__, __func__, __LINE__, msg);
#define HydraThrow(ExceptionType, message)                                     \
  throw ExceptionType("File    : " + std::string(__FILE__) + "\n" +            \
                      "       Function: " + std::string(__func__) + "\n" +     \
                      "       Line    : " + std::to_string(__LINE__) + "\n" +  \
                      "       " + std::string(message));
