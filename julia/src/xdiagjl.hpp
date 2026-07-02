// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include <utility>

#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <xdiag/all.hpp>

// The algorithm result structs are plain aggregates that CxxWrap would treat as
// "mirrored" types (mapping them field-by-field to a Julia struct) and then
// refuse to add_type. We instead expose them as opaque wrapped types with field
// accessors, so mirroring must be turned off for them here (visible in every
// translation unit that add_type's them or receives them as a return value).
namespace jlcxx {
template <> struct IsMirroredType<xdiag::EigsLanczosResult> : std::false_type {};
template <>
struct IsMirroredType<xdiag::EigvalsLanczosResult> : std::false_type {};
template <>
struct IsMirroredType<xdiag::EvolveLanczosResult> : std::false_type {};
template <>
struct IsMirroredType<xdiag::EvolveLanczosInplaceResult> : std::false_type {};
template <>
struct IsMirroredType<xdiag::TimeEvolveExpokitResult> : std::false_type {};
template <>
struct IsMirroredType<xdiag::TimeEvolveExpokitInplaceResult> : std::false_type {
};
template <> struct IsMirroredType<xdiag::EigsLobpcgResult> : std::false_type {};
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
