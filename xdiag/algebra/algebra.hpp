// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <set>
#include <string>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/blocks/boson.hpp>
#include <xdiag/blocks/electron.hpp>
#include <xdiag/blocks/fermion.hpp>
#include <xdiag/blocks/spinhalf.hpp>
#include <xdiag/blocks/tj.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>

namespace xdiag::algebra {

// Algebra describes how to bring an OpSum into normal order for a physical
// system. A block is rewritten by two plain functions (no std::function in the
// rules themselves -- only simple if/else):
//
//   expand(op, algebra)   one step of expanding a compound operator into
//                         simpler / elementary ones; nullopt if `op` is
//                         elementary or kept (see kept_named).
//   simplify(mono, alg)   one same-site / sorting / matrix-conversion step
//                         bringing a product into normal order; nullopt if it
//                         is already normal.
//
// normal_order() applies `expand` to a fixed point (skipping a size-1 operator
// whose type is in kept_named, so it survives for its dedicated matrix kernel),
// then `simplify` to a fixed point. The two rule bodies live in
// xdiag/algebra/rewrite/<block>_rules.{hpp,cpp}.
struct Algebra {
  std::string name; // human-readable name for errors
  int64_t nsites;   // number of sites
  int64_t d;        // local Hilbert space dimension per site
  std::set<std::string> fermionic_types; // anticommute at different sites
  std::set<std::string> allowed_types;   // valid input types (checked first)
  std::set<std::string> kept_named;      // size-1 ops kept named (kernels)
  std::function<std::optional<OpSum>(Op const &, Algebra const &)> expand;
  std::function<std::optional<OpSum>(Monomial const &, Algebra const &)>
      simplify;
};

// Concrete algebras (assembled in algebra.cpp). The "symmetry" variants reduce
// every operator to elementary generators in a site-major normal order (for the
// permutation-symmetry analysis); the "implementation" variants keep the named
// operators that have a dedicated matrix kernel and use the creation-major /
// matrix normal order the kernels expect.
Algebra spin_algebra(int64_t nsites);
Algebra matrix_algebra(int64_t nsites, int64_t d);
Algebra fermion_algebra(int64_t nsites);
Algebra electron_algebra(int64_t nsites);
Algebra tj_algebra(int64_t nsites);

Algebra fermion_implementation_algebra(int64_t nsites);
Algebra electron_implementation_algebra(int64_t nsites);
Algebra spinhalf_implementation_algebra(int64_t nsites);
// exchange_as_kernel = false drops Exchange from the kept (kernel) types so it
// expands to Cdag/C strings (used by the symmetric tJ block, which has no
// Exchange kernel).
Algebra tj_implementation_algebra(int64_t nsites, bool exchange_as_kernel = true);

// Per-block dispatch.
Algebra implementation_algebra(Boson const &block);
Algebra implementation_algebra(Spinhalf const &block);
Algebra implementation_algebra(Fermion const &block);
Algebra implementation_algebra(Electron const &block);
Algebra implementation_algebra(tJ const &block);
Algebra implementation_algebra(Block const &block);

Algebra symmetry_algebra(Boson const &block);
Algebra symmetry_algebra(Spinhalf const &block);
Algebra symmetry_algebra(Fermion const &block);
Algebra symmetry_algebra(Electron const &block);
Algebra symmetry_algebra(tJ const &block);
Algebra symmetry_algebra(Block const &block);

} // namespace xdiag::algebra
