// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "combine_matrix_ops.hpp"

#include <map>
#include <vector>

#include <xdiag/algebra/valid.hpp>
#include <xdiag/armadillo.hpp>
#include <xdiag/math/complex.hpp>
#include <xdiag/math/ipow.hpp>
#include <xdiag/utils/error.hpp>
#include <xdiag/utils/format.hpp>

namespace xdiag::operators {

// Embeds a single "Matrix" Op into the full dim_total × dim_total space, where
// dim_total = d^N and N = all_sites.size(). The op acts on a subset of size k
// of those sites; the remaining N-k sites must agree between row and column.
//
// Indexing: a global state index `r` decomposes as base-d digits, one per
// position in `all_sites`. The k digits at the op's positions form the local
// index into the op's matrix; the remaining N-k digits index the "untouched"
// subspace.
//
// Algorithm: iterate over (untouched, local_r, local_c) instead of (r, c).
// This is O(d^{N-k} * d^{2k}) work, equal to the number of non-zero entries,
// instead of the previous O(d^{2N}).
template <typename coeff_t>
static arma::Mat<coeff_t>
embed_op(Op const &op, std::vector<int64_t> const &all_sites,
         std::map<int64_t, int64_t> const &site_to_pos, int64_t d) {

  int64_t N = (int64_t)all_sites.size();
  int64_t dim_total = math::ipow(d, N);

  std::vector<int64_t> const &op_sites = op.sites();
  int64_t k = (int64_t)op_sites.size();
  int64_t dim_op = math::ipow(d, k);
  int64_t dim_other = math::ipow(d, N - k);

  // Positions of the op's sites within all_sites; mark which positions belong
  // to the op so we can enumerate the complementary "other" positions.
  std::vector<int64_t> op_pos(k);
  std::vector<bool> in_op(N, false);
  for (int64_t j = 0; j < k; ++j) {
    op_pos[j] = site_to_pos.at(op_sites[j]);
    in_op[op_pos[j]] = true;
  }
  std::vector<int64_t> other_pos;
  other_pos.reserve(N - k);
  for (int64_t p = 0; p < N; ++p) {
    if (!in_op[p]) {
      other_pos.push_back(p);
    }
  }

  // Precompute embeddings: local index → contribution to global index.
  std::vector<int64_t> op_stride(k);
  for (int64_t j = 0; j < k; ++j) {
    op_stride[j] = math::ipow(d, op_pos[j]);
  }
  std::vector<int64_t> other_stride(N - k);
  for (int64_t i = 0; i < N - k; ++i) {
    other_stride[i] = math::ipow(d, other_pos[i]);
  }

  std::vector<int64_t> local_to_global(dim_op, 0);
  for (int64_t local = 0; local < dim_op; ++local) {
    int64_t tmp = local;
    int64_t g = 0;
    for (int64_t j = 0; j < k; ++j) {
      g += (tmp % d) * op_stride[j];
      tmp /= d;
    }
    local_to_global[local] = g;
  }
  std::vector<int64_t> other_to_global(dim_other, 0);
  for (int64_t unt = 0; unt < dim_other; ++unt) {
    int64_t tmp = unt;
    int64_t g = 0;
    for (int64_t i = 0; i < N - k; ++i) {
      g += (tmp % d) * other_stride[i];
      tmp /= d;
    }
    other_to_global[unt] = g;
  }

  arma::Mat<coeff_t> const &M = op.matrix().as<arma::Mat<coeff_t>>();
  arma::Mat<coeff_t> full(dim_total, dim_total, arma::fill::zeros);
  for (int64_t unt = 0; unt < dim_other; ++unt) {
    int64_t base = other_to_global[unt];
    for (int64_t lr = 0; lr < dim_op; ++lr) {
      int64_t r = base + local_to_global[lr];
      for (int64_t lc = 0; lc < dim_op; ++lc) {
        int64_t c = base + local_to_global[lc];
        full(r, c) = M(lr, lc);
      }
    }
  }
  return full;
}

Op combine_matrix_ops(std::vector<Op> const &ops, int64_t d) try {
  if (ops.empty()) {
    XDIAG_THROW("combine_matrix_ops: ops vector must not be empty");
  }

  for (Op const &op : ops) {
    if (op.type() != "Matrix") {
      XDIAG_THROW(fmt::format(
          "combine_matrix_ops: all ops must be of type \"Matrix\", "
          "got \"{}\"",
          op.type()));
    }
    must_have_sites(op);
    must_have_matrix(op);
  }

  // Collect unique sites in first-appearance order
  std::vector<int64_t> all_sites;
  std::map<int64_t, int64_t> site_to_pos;
  for (Op const &op : ops) {
    for (int64_t s : op.sites()) {
      if (site_to_pos.find(s) == site_to_pos.end()) {
        site_to_pos[s] = (int64_t)all_sites.size();
        all_sites.push_back(s);
      }
    }
  }

  // Result is complex if any matrix is complex
  bool all_real = true;
  for (Op const &op : ops) {
    if (!op.matrix().isreal()) {
      all_real = false;
      break;
    }
  }

  // Build the product of embedded matrices
  if (all_real) {
    arma::mat combined = embed_op<double>(ops[0], all_sites, site_to_pos, d);
    for (std::size_t i = 1; i < ops.size(); ++i) {
      combined =
          combined * embed_op<double>(ops[i], all_sites, site_to_pos, d);
    }
    return Op("Matrix", all_sites, combined);
  } else {
    arma::cx_mat combined =
        embed_op<complex>(ops[0], all_sites, site_to_pos, d);
    for (std::size_t i = 1; i < ops.size(); ++i) {
      combined =
          combined * embed_op<complex>(ops[i], all_sites, site_to_pos, d);
    }
    return Op("Matrix", all_sites, combined);
  }
}
XDIAG_CATCH

} // namespace xdiag::operators
