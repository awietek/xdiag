// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstdint>

#include <xdiag/algebra/algebras/algebra.hpp>

namespace xdiag::algebra {

// Implementation algebra for the tJ block. Unlike tj_algebra (which reduces
// everything to elementary Cdagup/Cup/Cdagdn/Cdn for the symmetry analysis),
// this keeps the named operator types that have a dedicated matrix kernel (see
// matrices/blocks/tj) -- including the fast Exchange and the two diagonal
// spin-spin conventions SzSz and tJSzSz -- and only reduces the convenience
// types to them:
//   Hop      -> Hopup + Hopdn
//   Ntot     -> Nup + Ndn
//   Sz       -> 1/2 Nup - 1/2 Ndn
//   S+/S-    -> Cdagup Cdn / Cdagdn Cup ; Sx,Sy -> S+/S-
//   SdotS    -> SzSz   + Exchange   (plain spin coupling)
//   tJSdotS  -> tJSzSz + Exchange   (t-J coupling S.S - n n /4)
//   TotalN   -> sum_i Ntot{i}
// Exchange, SzSz, tJSzSz and the number operators are kept as size-1 named ops
// (routed to their kernels); they only expand to Cdag/C strings inside a product.
//
// exchange_as_kernel: the symmetric tJ block has no dedicated Exchange kernel --
// it reuses the electron Cdag/C string kernel -- so it passes false, which drops
// Exchange from the kernel types and lets it expand to (normal-ordered) Cdag/C
// strings like any other product. The non-symmetric block keeps the fast
// Exchange kernel (true).
Algebra tj_implementation_algebra(int64_t nsites,
                                  bool exchange_as_kernel = true);

} // namespace xdiag::algebra
