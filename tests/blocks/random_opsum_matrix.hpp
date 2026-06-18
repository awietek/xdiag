// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <algorithm>
#include <numeric>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <tests/catch.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/isapprox.hpp>
#include <xdiag/matrices/matrix.hpp>
#include <xdiag/operators/monomial.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/operators/types.hpp>
#include <xdiag/utils/logger.hpp>
#include <xdiag/utils/xdiag_show.hpp>

namespace xdiag::testcases {

// Build a random Op of the given type with random valid sites, honouring the
// op registry: site-free types carry no sites, and multi-site types draw
// distinct sites unless the type allows overlapping ones.
inline Op random_op(std::string const &type, int64_t nsites,
                    std::mt19937 &gen) {
  OpTypeInfo const &info = info_of_type(type);
  if (!info.site_required) { // Id, TotalN, HubbardU, ...
    return Op(type);
  }
  int64_t arity = info.nsites;
  std::uniform_int_distribution<int64_t> sitedist(0, nsites - 1);
  if (arity == 1) {
    return Op(type, sitedist(gen));
  }
  std::vector<int64_t> sites;
  if (info.allow_overlapping) {
    for (int64_t a = 0; a < arity; ++a) {
      sites.push_back(sitedist(gen));
    }
  } else {
    std::vector<int64_t> pool(nsites);
    std::iota(pool.begin(), pool.end(), int64_t(0));
    std::shuffle(pool.begin(), pool.end(), gen);
    sites.assign(pool.begin(), pool.begin() + arity);
  }
  return Op(type, sites);
}

// Shared randomized cross-check of the operator layer, applicable to any block:
// it builds a random OpSum (random complex coefficients, random monomial
// lengths, random sites) over every operator type the block supports, then
// verifies that the matrix produced by the full OpSum pipeline (expansion,
// same-site algebra, normal ordering, fermi signs, kernel dispatch) equals the
// matrix obtained by independently building each single-Op matrix and
// multiplying/adding them.
//
// The latter is the textbook definition of an operator polynomial, so equality
// validates the whole algebra machinery against first principles.
//
// The operator types are taken from the block's implementation algebra (so a
// new block is covered automatically), minus the explicit "Matrix" type (which
// needs an arma matrix) and any type whose site arity exceeds nsites.
//
// `block` MUST be on the full local Hilbert space (no conserved quantum
// numbers, no symmetry) so that every elementary Op is an endomorphism of the
// block and the matrix products are well defined.
template <typename block_t>
void test_random_opsum_matrix(block_t const &block, uint32_t seed,
                              int n_monomials = 40, int max_mono_size = 4) try {
  int64_t nsites = block.nsites();

  algebra::Algebra alg = algebra::implementation_algebra(block);
  std::vector<std::string> types;
  for (std::string const &t : alg.allowed_types) {
    OpTypeInfo const &info = info_of_type(t);
    if (info.matrix_required) { // "Matrix" needs an explicit arma matrix
      continue;
    }
    if (info.nsites != undefined && info.nsites > nsites) {
      continue; // not enough sites to place this operator
    }
    types.push_back(t);
  }

  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> cdist(-1.0, 1.0);
  std::uniform_int_distribution<int> sizedist(1, max_mono_size);
  std::uniform_int_distribution<std::size_t> typedist(0, types.size() - 1);

  OpSum opsum;
  std::vector<std::pair<complex, std::vector<Op>>> ref;
  for (int m = 0; m < n_monomials; ++m) {
    int size = sizedist(gen);
    std::vector<Op> ops;
    for (int s = 0; s < size; ++s) {
      ops.push_back(random_op(types[typedist(gen)], nsites, gen));
    }
    complex c(cdist(gen), cdist(gen));
    opsum += c * Monomial(ops);
    ref.emplace_back(c, ops);
  }

  // A: the full OpSum pipeline.
  arma::cx_mat A = matrixC(opsum, block);

  // B: independently multiply/add the individual single-Op matrices.
  int64_t d = A.n_rows;
  arma::cx_mat B(d, d, arma::fill::zeros);
  for (auto const &[c, ops] : ref) {
    arma::cx_mat prod = matrixC(ops.front(), block);
    for (std::size_t k = 1; k < ops.size(); ++k) {
      prod = prod * matrixC(ops[k], block);
    }
    B += c * prod;
  }

  REQUIRE(isapprox(A, B));
} catch (xdiag::Error const &e) {
  error_trace(e);
}

} // namespace xdiag::testcases
