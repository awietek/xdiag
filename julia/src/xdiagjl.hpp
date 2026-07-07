// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <cstdint>
#include <utility>

#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <xdiag/all.hpp>

using namespace xdiag;

// Register Mirrored types
namespace jlcxx {
template <> struct IsMirroredType<EigsLanczosResult> : std::false_type {};
template <> struct IsMirroredType<EigvalsLanczosResult> : std::false_type {};
template <> struct IsMirroredType<EigsLobpcgResult> : std::false_type {};
template <> struct IsMirroredType<EvolveLanczosResult> : std::false_type {};
template <>
struct IsMirroredType<EvolveLanczosInplaceResult> : std::false_type {};
template <> struct IsMirroredType<TimeEvolveExpokitResult> : std::false_type {};
template <>
struct IsMirroredType<TimeEvolveExpokitInplaceResult> : std::false_type {};
template <>
struct IsMirroredType<COOMatrix<int32_t, double>> : std::false_type {};
template <>
struct IsMirroredType<COOMatrix<int32_t, complex>> : std::false_type {};
template <>
struct IsMirroredType<COOMatrix<int64_t, double>> : std::false_type {};
template <>
struct IsMirroredType<COOMatrix<int64_t, complex>> : std::false_type {};
template <>
struct IsMirroredType<CSRMatrix<int32_t, double>> : std::false_type {};
template <>
struct IsMirroredType<CSRMatrix<int32_t, complex>> : std::false_type {};
template <>
struct IsMirroredType<CSRMatrix<int64_t, double>> : std::false_type {};
template <>
struct IsMirroredType<CSRMatrix<int64_t, complex>> : std::false_type {};
template <>
struct IsMirroredType<CSCMatrix<int32_t, double>> : std::false_type {};
template <>
struct IsMirroredType<CSCMatrix<int32_t, complex>> : std::false_type {};
template <>
struct IsMirroredType<CSCMatrix<int64_t, double>> : std::false_type {};
template <>
struct IsMirroredType<CSCMatrix<int64_t, complex>> : std::false_type {};
} // namespace jlcxx
#define JULIA_XDIAG_CALL_VOID(CMD)                                             \
  try {                                                                        \
    CMD;                                                                       \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_RETURN(CMD)                                           \
  try {                                                                        \
    return CMD;                                                                \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_RETURN_MOVE(CMD)                                      \
  try {                                                                        \
    return std::move(CMD);                                                     \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_ASSIGN(LVALUE, CMD)                                   \
  try {                                                                        \
    LVALUE = CMD;                                                              \
  } catch (xdiag::Error const &e) {                                            \
    xdiag::error_trace(e);                                                     \
    throw(std::runtime_error("Error occurred in XDiag C++ core library"));     \
  }
