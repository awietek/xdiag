// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

// Collector for the per-block dispatch_basis overloads. The variant-dispatch
// sites (kernels/apply.cpp, matrix.cpp, sparse/*.cpp) need every block's
// dispatch_basis overload visible at once (one body is instantiated for each
// alternative of the Block variant), so they include this single header instead
// of the individual ones. A new block adds its dispatch_basis.hpp here.
#include <xdiag/kernels/blocks/boson/dispatch_basis.hpp>
#include <xdiag/kernels/blocks/electron/dispatch_basis.hpp>
#include <xdiag/kernels/blocks/fermion/dispatch_basis.hpp>
#include <xdiag/kernels/blocks/spinhalf/dispatch_basis.hpp>
#include <xdiag/kernels/blocks/tj/dispatch_basis.hpp>
