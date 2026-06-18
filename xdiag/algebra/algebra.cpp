// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "algebra.hpp"

#include <variant>

#include <xdiag/algebra/rewrite/electron_rules.hpp>
#include <xdiag/algebra/rewrite/fermion_rules.hpp>
#include <xdiag/algebra/rewrite/matrix_rules.hpp>
#include <xdiag/algebra/rewrite/spin_rules.hpp>
#include <xdiag/algebra/rewrite/tj_rules.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/variants.hpp>

namespace xdiag::algebra {

// ---------------------------------------------------------------------------
// Concrete algebras. Each one only fills in the data of the Algebra struct: the
// allowed input types, the operators kept named for a matrix kernel, and the
// two plain rewrite functions (defined in rewrite/<block>_rules.cpp).
// ---------------------------------------------------------------------------

Algebra spin_algebra(int64_t nsites) {
  Algebra a;
  a.name = "spin-1/2";
  a.nsites = nsites;
  a.d = 2;
  a.allowed_types = {"Exchange",        "ExchangeAsym", "Id",   "S+",
                     "S-",              "ScalarChirality", "SdotS", "Sx",
                     "Sy",              "Sz",           "SzSz", "TotalSz"};
  a.expand = spin_expand;
  a.simplify = spin_simplify;
  return a;
}

Algebra matrix_algebra(int64_t nsites, int64_t d) {
  Algebra a;
  a.name = "matrix";
  a.nsites = nsites;
  a.d = d;
  a.allowed_types = {"A",      "Adag",  "Exchange", "ExchangeAsym",
                     "Hop",    "HopAsym", "HubbardU", "Id",
                     "Matrix", "N",     "S+",       "S-",
                     "ScalarChirality", "SdotS", "Sx", "Sy",
                     "Sz",     "SzSz",  "TotalN"};
  // Nothing is kept named: every operator is converted to an explicit Matrix.
  a.expand = boson_expand;
  a.simplify = matrix_simplify;
  return a;
}

Algebra fermion_algebra(int64_t nsites) {
  Algebra a;
  a.name = "fermion";
  a.nsites = nsites;
  a.d = 2;
  a.fermionic_types = {"C", "Cdag"};
  a.allowed_types = {"C", "Cdag", "Hop", "HopAsym", "Id", "N", "NN", "TotalN"};
  // Empty kept set: everything expands to elementary C / Cdag.
  a.expand = fermion_expand;
  a.simplify = fermion_simplify;
  return a;
}

Algebra fermion_implementation_algebra(int64_t nsites) {
  Algebra a = fermion_algebra(nsites);
  a.name = "fermion_implementation";
  a.kept_named = {"Hop", "HopAsym", "N", "NN", "TotalN"};
  return a;
}

Algebra electron_algebra(int64_t nsites) {
  Algebra a;
  a.name = "electron";
  a.nsites = nsites;
  a.d = 4;
  a.fermionic_types = {"Cdagup", "Cup", "Cdagdn", "Cdn"};
  a.allowed_types = {"Cdagdn",   "Cdagup", "Cdn",      "Cup",
                     "Exchange", "ExchangeAsym", "Hop",  "HopAsym",
                     "Hopdn",    "HopdnAsym", "Hopup",   "HopupAsym",
                     "HubbardU", "Id",     "Ndn",      "NdnNdn",
                     "NdnNup",   "Ntot",   "NtotNtot", "Nup",
                     "NupNdn",   "NupNup", "Nupdn",    "NupdnNupdn",
                     "S+",       "S-",     "SdotS",    "Sx",
                     "Sy",       "Sz",     "SzSz",     "TotalN",
                     "TotalNup", "TotalNdn", "TotalSz"};
  a.expand = electron_expand;
  a.simplify = electron_simplify_symmetry;
  return a;
}

Algebra electron_implementation_algebra(int64_t nsites) {
  Algebra a = electron_algebra(nsites);
  a.name = "electron_implementation";
  a.kept_named = {"Cdagup", "Cup",        "Cdagdn",   "Cdn",    "Hopup",
                  "Hopdn",  "HopupAsym",  "HopdnAsym", "Nup",   "Ndn",
                  "Nupdn",  "NupdnNupdn", "NtotNtot",  "NupNdn", "NupNup",
                  "NdnNdn", "NdnNup",     "SzSz",      "HubbardU"};
  a.simplify = electron_simplify; // creation-major
  return a;
}

Algebra tj_algebra(int64_t nsites) {
  Algebra a;
  a.name = "tJ";
  a.nsites = nsites;
  a.d = 3;
  a.fermionic_types = {"Cdagup", "Cup", "Cdagdn", "Cdn"};
  a.allowed_types = {"Cdagdn",   "Cdagup", "Cdn",      "Cup",
                     "Exchange", "ExchangeAsym", "Hop",  "HopAsym",
                     "Hopdn",    "HopdnAsym", "Hopup",   "HopupAsym",
                     "Id",       "Ndn",    "NdnNdn",   "NdnNup",
                     "Ntot",     "NtotNtot", "Nup",    "NupNdn",
                     "NupNup",   "S+",     "S-",       "SdotS",
                     "SzSz",     "Sx",     "Sy",       "Sz",
                     "TotalN",   "TotalNup", "TotalNdn", "TotalSz",
                     "tJSdotS",  "tJSzSz"};
  a.expand = tj_expand;
  a.simplify = tj_simplify_symmetry;
  return a;
}

Algebra tj_implementation_algebra(int64_t nsites, bool exchange_as_kernel) {
  Algebra a = tj_algebra(nsites);
  a.name = "tj_implementation";
  a.kept_named = {"Cdagup", "Cup",    "Cdagdn",   "Cdn",    "Hopup",
                  "Hopdn",  "Nup",    "Ndn",      "NupNup", "NdnNdn",
                  "NupNdn", "NdnNup", "NtotNtot", "SzSz",   "tJSzSz"};
  if (exchange_as_kernel) {
    a.kept_named.insert("Exchange");
  }
  a.simplify = tj_simplify; // creation-major, projected
  return a;
}

Algebra spinhalf_implementation_algebra(int64_t nsites) {
  Algebra a;
  a.name = "spinhalf_implementation";
  a.nsites = nsites;
  a.d = 2;
  a.allowed_types = {"Exchange", "ExchangeAsym", "Id", "Matrix", "S+", "S-",
                     "ScalarChirality", "SdotS", "Sx", "Sy", "Sz", "SzSz",
                     "TotalSz"};
  // Named operators with a dedicated spin-1/2 kernel; the rest are converted to
  // local Matrix operators.
  a.kept_named = {"Exchange", "ExchangeAsym", "S+",  "S-",
                  "ScalarChirality", "Sz",    "SzSz"};
  a.expand = spinhalf_expand;
  a.simplify = matrix_simplify;
  return a;
}

// ---------------------------------------------------------------------------
// Per-block dispatch.
// ---------------------------------------------------------------------------

Algebra implementation_algebra(Spinhalf const &block) try {
  return spinhalf_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Boson const &block) try {
  return matrix_algebra(block.nsites(), block.d());
}
XDIAG_CATCH

Algebra implementation_algebra(Fermion const &block) try {
  return fermion_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Electron const &block) try {
  return electron_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(tJ const &block) try {
  return tj_implementation_algebra(block.nsites());
}
XDIAG_CATCH

Algebra implementation_algebra(Block const &block) try {
  return std::visit([](auto const &b) { return implementation_algebra(b); },
                    block);
}
XDIAG_CATCH

Algebra symmetry_algebra(Spinhalf const &block) try {
  return matrix_algebra(block.nsites(), 2);
}
XDIAG_CATCH

Algebra symmetry_algebra(Boson const &block) try {
  return matrix_algebra(block.nsites(), block.d());
}
XDIAG_CATCH

Algebra symmetry_algebra(Fermion const &block) try {
  return fermion_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(Electron const &block) try {
  return electron_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(tJ const &block) try {
  return tj_algebra(block.nsites());
}
XDIAG_CATCH

Algebra symmetry_algebra(Block const &block) try {
  return std::visit([](auto const &b) { return symmetry_algebra(b); }, block);
}
XDIAG_CATCH

} // namespace xdiag::algebra
