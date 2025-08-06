// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <exception>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <xdiag/config.hpp>
#include <xdiag/utils/xdiag_api.hpp>

namespace xdiag {

class Error // : public std::exception
{
public:
  explicit Error(std::string message, std::string original_message = "");
  Error(Error const &error, std::string message);
  // virtual XDIAG_API ~Error();
  std::vector<std::string> const &messages() const;
  std::string original_message() const;
  XDIAG_API const char *what() const noexcept;

private:
  std::string original_message_;
  std::vector<std::string> messages_;
};

XDIAG_API void throw_error(std::string message, const char *file,
                           const char *func, int line);
XDIAG_API void rethrow_error(Error const &error, const char *file,
                             const char *func, int line);
XDIAG_API void error_trace(Error const &error);

std::string cut_file(const char *file);

// common error checks
void check_dimension_works_with_blas_int_size(int64_t dim);
template <typename bit_t> void check_nsites_work_with_bits(int64_t nsites);

} // namespace xdiag

// comment: throw 0; necessary to let compiler know there need not be a return
// value
#define XDIAG_THROW(message)                                                   \
  xdiag::throw_error(message, __FILE__, __func__, __LINE__);                   \
  throw 0;

#define XDIAG_RETHROW(error)                                                   \
  xdiag::rethrow_error(error, __FILE__, __func__, __LINE__);                   \
  throw 0;

#define XDIAG_CATCH                                                            \
  catch (Error const &e) {                                                     \
    xdiag::rethrow_error(e, __FILE__, __func__, __LINE__);                     \
    throw 0;                                                                   \
  }
#define XDIAG_TRY_CATCH(statement)                                             \
  try {                                                                        \
    statement;                                                                 \
  } catch (Error const &e) {                                                   \
    xdiag::rethrow_error(e, __FILE__, __func__, __LINE__);                     \
    throw 0;                                                                   \
  }
