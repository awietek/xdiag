// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

// Dense cross-checks of symmetry-adapted blocks against the plain block: the
// spectrum of the full block must equal the union of the spectra of its
// symmetry sectors. Any missing or double-counted representative, or a wrong
// phase in a symmetric kernel, shows up here.

#include <algorithm>
#include <functional>
#include <vector>

#include <tests/catch.hpp>
#include <tests/is_approx_hermitian.hpp>

#include <xdiag/blocks/blocks.hpp>
#include <xdiag/kernels/matrix.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/symmetries/cyclic_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::testcases {

// Dense spectrum of ops on a block via the complex matrix (works for every
// momentum sector), with a hermiticity check. For a real block and real ops the
// real kernel must give the same spectrum.
inline arma::vec dense_spectrum(OpSum const &ops, Block const &block) {
  arma::cx_mat H = matrixC(ops, block);
  REQUIRE(is_approx_hermitian(H, 1e-8));
  arma::vec e;
  arma::eig_sym(e, H);
  if (isreal(block) && isreal(ops)) {
    arma::mat H_real = matrix(ops, block);
    arma::vec e_real;
    arma::eig_sym(e_real, H_real);
    REQUIRE(isapprox(e, e_real));
  }
  return e;
}

// The spectrum of `full` must equal the union of the spectra of `parts`, each
// repeated `mults` times (irrep dimension). Empty parts are skipped.
inline void spectrum_union(OpSum const &ops, Block const &full,
                           std::vector<Block> const &parts,
                           std::vector<int64_t> const &mults = {}) {
  arma::vec e_full = dense_spectrum(ops, full);

  std::vector<double> e_union;
  int64_t dimsum = 0;
  for (int64_t k = 0; k < (int64_t)parts.size(); ++k) {
    int64_t mult = mults.empty() ? 1 : mults[k];
    dimsum += dim(parts[k]) * mult;
    if (dim(parts[k]) == 0) {
      continue;
    }
    arma::vec e_k = dense_spectrum(ops, parts[k]);
    for (double e : e_k) {
      for (int64_t m = 0; m < mult; ++m) {
        e_union.push_back(e);
      }
    }
  }
  REQUIRE(dimsum == dim(full));
  std::sort(e_union.begin(), e_union.end());
  REQUIRE(isapprox(e_full, arma::vec(e_union)));
}

// spectrum_union over all momentum sectors of the cyclic group on `nsites`
// sites; `sector(irrep)` returns the symmetry-adapted block for an irrep.
inline void translation_spectrum_union(
    OpSum const &ops_pbc, Block const &full, int64_t nsites,
    std::function<Block(Representation const &)> const &sector) {
  std::vector<Block> parts;
  for (int64_t k = 0; k < nsites; ++k) {
    parts.push_back(sector(cyclic_group_irrep(nsites, k)));
  }
  spectrum_union(ops_pbc, full, parts);
}

} // namespace xdiag::testcases
