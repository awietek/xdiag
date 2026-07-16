// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <string>

#include <xdiag/basis/basis.hpp>
#include <xdiag/utils/error.hpp>

namespace xdiag::kernels {

// Resolves two type-erased Basis objects to one concrete basis type from the
// list Bs... and invokes f(basis_in, basis_out) with the concrete references.
//
// The concrete type is found by matching the runtime basis type-id against each
// Bs::static_type() in a compile-time fold: N integer comparisons, short-circuited
// at the first match, no allocation and no std::function. Order the hot bases
// first in the list if it matters. Throws if the two bases differ in type or if
// none of Bs... matches (the latter means the caller listed a basis here but did
// not instantiate the kernel for it, or vice versa).
template <typename... Bs, typename func_t>
void dispatch_basis_types(basis::Basis const &a, basis::Basis const &b,
                          func_t &&f) try {
  if (a.type() != b.type()) {
    XDIAG_THROW("Type mismatch for Basis: \n  1 -> " + std::string(a.name()) +
                "\n  2 -> " + std::string(b.name()));
  }
  std::size_t t = a.type();
  // Unary right fold over `||` (the `( pattern || ... )` syntax). `Bs` is a
  // parameter pack; the `...` repeats the pattern once per type in Bs and joins
  // the terms with `||`, so for Bs = {B1, B2, ...} this expands to:
  //
  //   matched = (t == B1::static_type() ? (f(cast<B1>(a), cast<B1>(b)), true)
  //                                      : false)
  //          || (t == B2::static_type() ? (f(cast<B2>(a), cast<B2>(b)), true)
  //                                      : false)
  //          || ... ;
  //
  // i.e. a compiler-generated if/else-if chain, one branch per listed type. Per
  // term: if the runtime id `t` matches this concrete type's static_type(), the
  // comma operator `(f(...), true)` runs f (the static_cast is safe because the
  // id matched) and yields true; otherwise the term is false. `||` short-circuits,
  // so f runs exactly once at the first match and the rest are skipped. `matched`
  // is false only if no type in Bs matched.
  bool matched = ((t == Bs::static_type()
                       ? (f(static_cast<Bs const &>(a),
                           static_cast<Bs const &>(b)),
                          true)
                       : false) ||
                  ...);
  if (!matched) {
    XDIAG_THROW("No handler registered for type: " + std::string(a.name()));
  }
}
XDIAG_CATCH

} // namespace xdiag::kernels
