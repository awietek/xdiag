#pragma once

#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <xdiag/config.hpp>

namespace xdiag {

class Error : public std::exception {
public:
  explicit Error(std::string message, std::string original_message = "");
  Error(Error const &error, std::string message);
  std::vector<std::string> const &messages() const;
  std::string original_message() const;
  const char *what() const noexcept;

private:
  std::string original_message_;
  std::vector<std::string> messages_;
};

void throw_error(std::string message, const char *file, const char *func,
                 int line);
void rethrow_error(Error const &error, const char *file, const char *func,
                   int line);
void error_trace(Error const &error);

std::string cut_file(const char *file);

#define XDIAG_THROW(message)                                                   \
  xdiag::throw_error(message, __FILE__, __func__, __LINE__);

#define XDIAG_RETHROW(error)                                                   \
  xdiag::rethrow_error(error, __FILE__, __func__, __LINE__);

template <class... Args>
inline void rethrow2(const char *file, const char *func, int line,
                     Args &&...args) {
  std::cout << "RT1\n";

  // build an error message
  std::ostringstream ss;
  ss << "File    : " << cut_file(file) << "\n";
  ss << "       Func    : " << func << " (line " << line << ")\n";
  ss << "       ";
  auto sep = " : ";
  (void)sep;
  std::cout << "RT2\n";
  using expand = int[];
  void(expand{0, ((ss << args), sep = ", ", 0)...});
  // figure out what kind of exception is active
  std::cout << "RT3\n";
  try {
    std::cout << "RTrt\n";
    std::rethrow_exception(std::current_exception());
  } // catch (const std::runtime_error &e) {
  //   std::cout << "RTru\n";
  //   std::throw_with_nested(std::runtime_error(ss.str()));
  // } catch (const std::invalid_argument &e) {
  //   std::cout << "RTia\n";
  //   std::throw_with_nested(std::invalid_argument(ss.str()));
  // } catch (const std::logic_error &e) {
  //   std::cout << "RTle\n";
  //   std::throw_with_nested(std::logic_error(ss.str()));
  // }  catch (const std::nested_exception &e) {
  //   std::cout << "RTnest\n";
  //   std::throw_with_nested(std::logic_error(ss.str()));
  // }
  // etc - default to a runtime_error
  catch (Error const &e) {
    std::cout << "RTerr\n";
    std::throw_with_nested(std::runtime_error(ss.str()));
  }

  catch (...) {
    std::cout << "RT4df\n";
    // std::throw_with_nested(std::runtime_error(ss.str()));
    // throw std::runtime_error(ss.str());
    std::throw_with_nested(std::runtime_error(ss.str()));
    // std::throw_with_nested(ss.str());
  }
  std::cout << "RT5\n";
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

} // namespace xdiag

#define XDiagRethrow(msg) xdiag::rethrow2(__FILE__, __func__, __LINE__, msg);
#define XDiagThrow(ExceptionType, message)                                     \
  throw ExceptionType("File    : " + cut_file(__FILE__) + "\n" +               \
                      "       In      : " + std::string(__func__) +            \
                      " (line " + std::to_string(__LINE__) + ")\n" +           \
                      "       " + std::string(message));
