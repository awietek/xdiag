// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

// Visibility control
#ifdef XDIAG_STATIC_DEFINE
#define XDIAG_API
#else
#include <xdiag/export.hpp>

#if defined(XDIAG_EXPORT)
#define XDIAG_API XDIAG_EXPORT
#elif defined(XDIAG_DISTRIBUTED_EXPORT)
#define XDIAG_API XDIAG_DISTRIBUTED_EXPORT
#elif defined(XDIAGJL_EXPORT)
#define XDIAG_API XDIAGJL_EXPORT
#else
#error "Cannot determine library name for export symbols"
#endif

#endif
