// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/matrices/blocks/boson/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/boson/matrix_policy.hpp>
#include <xdiag/matrices/blocks/fermion/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/fermion/matrix_policy.hpp>
#include <xdiag/matrices/blocks/spinhalf/dispatch_basis.hpp>
#include <xdiag/matrices/blocks/spinhalf/matrix_policy.hpp>

namespace xdiag::matrices {

// Registry mapping a concrete Block type to the pieces the generic matrix /
// apply / sparse kernels need:
//
//   policy   — the block's MatrixPolicy (wraps its matrix_generic kernel).
//   dispatch — resolves the block's shared_ptr<Basis> to concrete
//              BasisOnTheFly<...> / BasisSymmetric<...> objects and invokes
//              f(basis_in, basis_out).
//
// This is the single place that knows how to drive a block's numerical kernel.
// The generic builders are templates over block_t that route through here; the
// public Block-variant entry points unwrap the variant exactly once (via
// visit_same_type) and call those templates.
//
// Registering a NEW block type:
//   1. add a kernel_traits specialization below,
//   2. add it to the Block variant (blocks.hpp), and
//   3. add a kernel instantiation group for it (matrices/kernels.cpp).
// The primary template is intentionally undefined: a missing specialization is
// a compile error at the dispatch site (not a silent fallback that recurses).
template <typename block_t> struct kernel_traits {
  static_assert(
      sizeof(block_t) == 0,
      "kernel_traits is not specialized for this block type. Register "
      "the block by adding a kernel_traits specialization here, a "
      "Block variant alternative, and a kernel instantiation group.");
};

template <> struct kernel_traits<Spinhalf> {
  using policy = spinhalf::MatrixPolicy;
  template <typename func_t>
  static void dispatch(Spinhalf const &block_in, Spinhalf const &block_out,
                       func_t &&f) {
    spinhalf::dispatch_basis(block_in, block_out, std::forward<func_t>(f));
  }
};

template <> struct kernel_traits<Boson> {
  using policy = boson::MatrixPolicy;
  template <typename func_t>
  static void dispatch(Boson const &block_in, Boson const &block_out,
                       func_t &&f) {
    boson::dispatch_basis(block_in, block_out, std::forward<func_t>(f));
  }
};

template <> struct kernel_traits<Fermion> {
  using policy = fermion::MatrixPolicy;
  template <typename func_t>
  static void dispatch(Fermion const &block_in, Fermion const &block_out,
                       func_t &&f) {
    fermion::dispatch_basis(block_in, block_out, std::forward<func_t>(f));
  }
};

} // namespace xdiag::matrices
