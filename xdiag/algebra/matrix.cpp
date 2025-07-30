// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "matrix.hpp"

#include <xdiag/algebra/apply_dispatch.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {

template <typename coeff_t, typename op_t, class block_t>
static arma::Mat<coeff_t> matrix_gen(op_t const &op, block_t const &block_in,
                                     block_t const &block_out) try {
  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix(mat.memptr(), OpSum(op), block_in, block_out);
  return mat;
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, typename op_t, typename block_t>
static arma::Mat<coeff_t> matrix_gen(op_t const &op,
                                     block_t const &blocki) try {
  auto blockr = block(OpSum(op), blocki);
  return matrix_gen<coeff_t>(op, blocki, blockr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, typename op_t>
static arma::Mat<coeff_t> matrix_gen_variant(op_t const &op,
                                             Block const &block_in,
                                             Block const &block_out) try {
  return std::visit(
      overload{
          [&](Spinhalf const &b1, Spinhalf const &b2) {
            return matrix_gen<coeff_t>(op, b1, b2);
          },
          [&](tJ const &b1, tJ const &b2) {
            return matrix_gen<coeff_t>(op, b1, b2);
          },
          [&](Electron const &b1, Electron const &b2) {
            return matrix_gen<coeff_t>(op, b1, b2);
          },
#ifdef XDIAG_USE_MPI
          [&](SpinhalfDistributed const &, SpinhalfDistributed const &) {
            XDIAG_THROW("Matrix creation not implemented for "
                        "SpinhalfDistributed blocks");
            return arma::Mat<coeff_t>();
          },
          [&](tJDistributed const &, tJDistributed const &) {
            XDIAG_THROW(
                "Matrix creation not implemented for tJDistributed blocks");
            return arma::Mat<coeff_t>();
          },
          [&](ElectronDistributed const &, ElectronDistributed const &) {
            XDIAG_THROW("Matrix creation not implemented for "
                        "ElectronDistributed blocks");
            return arma::Mat<coeff_t>();
          },
#endif
          [&](auto &&, auto &&) {
            XDIAG_THROW(fmt::format("Invalid combination of Block types"));
            return arma::Mat<coeff_t>();
          }},
      block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename coeff_t, typename op_t>
static arma::Mat<coeff_t> matrix_gen_variant(op_t const &op,
                                             Block const &blocki) try {
  auto blocko = block(OpSum(op), blocki);
  return matrix_gen_variant<coeff_t>(op, blocki, blocko);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

arma::mat matrix(Op const &op, Block const &block) try {
  return std::visit(
      [&](auto const &b) { return matrix_gen_variant<double>(op, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::mat matrix(OpSum const &ops, Block const &block) try {
  return std::visit(
      [&](auto const &b) { return matrix_gen_variant<double>(ops, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(Op const &op, Block const &block) try {
  return std::visit(
      [&](auto const &b) { return matrix_gen_variant<complex>(op, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(OpSum const &ops, Block const &block) try {
  return std::visit(
      [&](auto const &b) { return matrix_gen_variant<complex>(ops, b); },
      block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

arma::mat matrix(Op const &op, Block const &block_in,
                 Block const &block_out) try {
  return matrix_gen_variant<double>(op, block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::mat matrix(OpSum const &ops, Block const &block_in,
                 Block const &block_out) try {
  return matrix_gen_variant<double>(ops, block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(Op const &op, Block const &block_in,
                     Block const &block_out) try {
  return matrix_gen_variant<complex>(op, block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(OpSum const &ops, Block const &block_in,
                     Block const &block_out) try {
  return matrix_gen_variant<complex>(ops, block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// developer methods
template <typename coeff_t, class block_t>
void matrix(coeff_t *mat, OpSum const &ops, block_t const &block_in,
            block_t const &block_out) try {

  // Check if ops and blocks are compatible
  if (!blocks_match(ops, block_in, block_out)) {
    XDIAG_THROW("Cannot matrix on Blocks. The resulting Block is not in "
                "the correct symmetry sector. Please check the quantum numbers "
                "of the output block.");
  }

  // Check if real matrix can be created
  if constexpr (isreal<coeff_t>()) {
    if (!isreal(ops)) {
      XDIAG_THROW(
          "Cannot create a real matrix from an Op or OpSum which is complex. "
          "Please use the function \"matrixC\" instead.");
    }
    if (!isreal(block_in) || !isreal(block_out)) {
      XDIAG_THROW("Cannot create a real matrix when a block is complex. "
                  "Please use the function \"matrixC\" instead.")
    }
  }

  int64_t nsites = block_in.nsites();
  check_valid(ops, nsites);
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);
  OpSum opsc = operators::compile<block_t>(ops);

  // create fill method to add to the matrix
  auto fill = [&](int64_t idx_in, int64_t idx_out, coeff_t val) {
    mat[idx_out + idx_in * m] += val;
  };

#ifdef _OPENMP
  omp_set_schedule(omp_sched_guided, 0);
#endif
  
  algebra::apply_dispatch<coeff_t>(opsc, block_in, block_out, fill);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template XDIAG_API void matrix(double *, OpSum const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(double *, OpSum const &, tJ const &, tJ const &);
template XDIAG_API void matrix(double *, OpSum const &, Electron const &,
                               Electron const &);
template XDIAG_API void matrix(complex *, OpSum const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(complex *, OpSum const &, tJ const &,
                               tJ const &);
template XDIAG_API void matrix(complex *, OpSum const &, Electron const &,
                               Electron const &);

} // namespace xdiag
