// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

// Branch-prediction hints.
// [[likely]] / [[unlikely]] are C++20; use __builtin_expect on C++17.
//
// Compiler coverage:
//   GCC            : defines __GNUC__
//   Clang          : defines __clang__
//   Intel icpx     : defines __clang__ and __INTEL_LLVM_COMPILER
//   Intel icpc     : defines __INTEL_COMPILER (and __GNUC__ for compat)
//   MSVC / unknown : falls back to no-op
#if defined(__GNUC__) || defined(__clang__) || \
    defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define XDIAG_LIKELY(x)   (__builtin_expect(!!(x), 1))
#define XDIAG_UNLIKELY(x) (__builtin_expect(!!(x), 0))
#else
#define XDIAG_LIKELY(x)   (x)
#define XDIAG_UNLIKELY(x) (x)
#endif
