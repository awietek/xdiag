#include "matrix.hpp"

#include <xdiag/algebra/fill.hpp>
#include <xdiag/basis/electron/apply/dispatch_matrix.hpp>
#include <xdiag/basis/spinhalf/apply/dispatch_matrix.hpp>
#include <xdiag/basis/tj/apply/dispatch_matrix.hpp>
#include <xdiag/operators/logic/block.hpp>
#include <xdiag/operators/logic/compilation.hpp>
#include <xdiag/operators/logic/real.hpp>
#include <xdiag/operators/logic/valid.hpp>

namespace xdiag {

template <typename coeff_t, typename op_t, class block_t>
static arma::Mat<coeff_t> matrix_gen(op_t const &op, block_t const &block_in,
                                     block_t const &block_out) try {
  if (!blocks_match(OpSum(op), block_in, block_out)) {
    XDIAG_THROW("Cannot matrix on Blocks. The resulting Block is not in "
                "the correct symmetry sector. Please check the quantum numbers "
                "of the output block.");
  }
  int64_t m = block_out.size();
  int64_t n = block_in.dim();
  arma::Mat<coeff_t> mat(m, n, arma::fill::zeros);
  matrix(mat.memptr(), op, block_in, block_out);
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

  // Check if real matrix can be created
  if constexpr (isreal<coeff_t>()) {
    if (!isreal(op)) {
      XDIAG_THROW(
          "Cannot create a real matrix from an Op or OpSum which is complex. "
          "Please use the function \"matrixC\" instead.");
    }
    if (!isreal(block_in) || !isreal(block_out)) {
      XDIAG_THROW("Cannot create a real matrix when a block is complex. "
                  "Please use the function \"matrixC\" instead.")
    }
  }

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
          [](SpinhalfDistributed const &, SpinhalfDistributed const &) {
            XDIAG_THROW("Matrix creation not implemented for "
                        "SpinhalfDistributed blocks");
            return arma::Mat<coeff_t>();
          },
          [](tJ const &, tJ const &) {
            XDIAG_THROW(
                "Matrix creation not implemented for tJDistributed blocks");
            return arma::Mat<coeff_t>();
          },
#endif
          [](auto const &, auto const &) {
            XDIAG_THROW(fmt::format("Invalid combination of Block types"));
            return arma::Mat<coeff_t>();
          }},
      block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

arma::mat matrix(Op const &op, Block const &block) try {
  return std::visit([&](auto const &b) { return matrix_gen<double>(op, b); },
                    block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::mat matrix(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return matrix_gen<double>(ops, b); },
                    block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(Op const &op, Block const &block) try {
  return std::visit([&](auto const &b) { return matrix_gen<complex>(op, b); },
                    block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
arma::cx_mat matrixC(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return matrix_gen<complex>(ops, b); },
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

template <typename coeff_t, class block_t>
static void compile_and_dispatch(coeff_t *mat, int64_t m, OpSum const &ops,
                                 block_t const &block_in,
                                 block_t const &block_out);

template <>
void compile_and_dispatch(double *mat, int64_t m, OpSum const &ops,
                          Spinhalf const &block_in,
                          Spinhalf const &block_out) try {
  OpSum opsc = operators::compile_spinhalf(ops);
  basis::spinhalf::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void compile_and_dispatch(complex *mat, int64_t m, OpSum const &ops,
                          Spinhalf const &block_in,
                          Spinhalf const &block_out) try {
  OpSum opsc = operators::compile_spinhalf(ops);
  basis::spinhalf::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void compile_and_dispatch(double *mat, int64_t m, OpSum const &ops,
                          tJ const &block_in, tJ const &block_out) try {
  OpSum opsc = operators::compile_tj(ops);
  basis::tj::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void compile_and_dispatch(complex *mat, int64_t m, OpSum const &ops,
                          tJ const &block_in, tJ const &block_out) try {
  OpSum opsc = operators::compile_tj(ops);
  basis::tj::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
void compile_and_dispatch(double *mat, int64_t m, OpSum const &ops,
                          Electron const &block_in,
                          Electron const &block_out) try {
  OpSum opsc = operators::compile_electron(ops);
  basis::electron::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template <>
void compile_and_dispatch(complex *mat, int64_t m, OpSum const &ops,
                          Electron const &block_in,
                          Electron const &block_out) try {
  OpSum opsc = operators::compile_electron(ops);
  basis::electron::dispatch_matrix(opsc, block_in, block_out, mat, m);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

// developer methods
template <typename coeff_t, class block_t>
void matrix(coeff_t *mat, OpSum const &ops, block_t const &block_in,
            block_t const &block_out) try {
  int64_t nsites = block_in.nsites();
  check_valid(ops, nsites);
  int64_t m = block_out.size();
  int64_t n = block_in.size();
  std::fill(mat, mat + m * n, 0);
  compile_and_dispatch(mat, m, ops, block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template XDIAG_API void matrix(double *, OpSum const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(complex *, OpSum const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(double *, OpSum const &, tJ const &, tJ const &);
template XDIAG_API void matrix(complex *, OpSum const &, tJ const &,
                               tJ const &);
template XDIAG_API void matrix(double *, OpSum const &, Electron const &,
                               Electron const &);
template XDIAG_API void matrix(complex *, OpSum const &, Electron const &,
                               Electron const &);

template <typename coeff_t, class block_t>
void matrix(coeff_t *mat, Op const &op, block_t const &block_in,
            block_t const &block_out) try {
  matrix(mat, OpSum(op), block_in, block_out);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template XDIAG_API void matrix(double *, Op const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(complex *, Op const &, Spinhalf const &,
                               Spinhalf const &);
template XDIAG_API void matrix(double *, Op const &, tJ const &, tJ const &);
template XDIAG_API void matrix(complex *, Op const &, tJ const &, tJ const &);
template XDIAG_API void matrix(double *, Op const &, Electron const &,
                               Electron const &);
template XDIAG_API void matrix(complex *, Op const &, Electron const &,
                               Electron const &);

template <typename coeff_t, class block_t>
void matrix(coeff_t *mat, OpSum const &ops, block_t const &blocki) try {
  auto blockr = block(ops, blocki);
  matrix(mat, ops, blocki, blockr);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template XDIAG_API void matrix(double *, OpSum const &, Spinhalf const &);
template XDIAG_API void matrix(complex *, OpSum const &, Spinhalf const &);
template XDIAG_API void matrix(double *, OpSum const &, tJ const &);
template XDIAG_API void matrix(complex *, OpSum const &, tJ const &);
template XDIAG_API void matrix(double *, OpSum const &, Electron const &);
template XDIAG_API void matrix(complex *, OpSum const &, Electron const &);

template <typename coeff_t, class block_t>
void matrix(coeff_t *mat, Op const &op, block_t const &blocki) try {
  matrix(mat, OpSum(op), blocki);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template XDIAG_API void matrix(double *, Op const &, Spinhalf const &);
template XDIAG_API void matrix(complex *, Op const &, Spinhalf const &);
template XDIAG_API void matrix(double *, Op const &, tJ const &);
template XDIAG_API void matrix(complex *, Op const &, tJ const &);
template XDIAG_API void matrix(double *, Op const &, Electron const &);
template XDIAG_API void matrix(complex *, Op const &, Electron const &);

} // namespace xdiag
