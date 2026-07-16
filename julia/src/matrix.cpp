// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written special -- copied verbatim by generate.sh.
//
// Dense matrix builder: fills a caller-allocated, Julia-owned column-major
// buffer of size size(out) x dim(in). The output block is inferred via
// block(ops, in) (so Julia never round-trips the Block variant). Templated on
// coeff_t, so instantiated here per block. Sparse builders live in sparse.cpp.

#include <julia/src/xdiagjl.hpp>

namespace xdiag::julia {

// OpSum: zero-copy pointer-fill into a Julia-owned buffer. Op / Monomial reuse
// this via matrix(OpSum(op), block) on the Julia side (a one-term OpSum).
template <typename block_t> static void register_matrix(jlcxx::Module &mod) {
  mod.method("fun_matrix", [](OpSum const &ops, block_t const &in, double *m) {
    JULIA_XDIAG_CALL_VOID((matrix<double>(ops, in, block(ops, in), m)));
  });
  mod.method(
      "fun_matrixC", [](OpSum const &ops, block_t const &in, complex *m) {
        JULIA_XDIAG_CALL_VOID((matrix<complex>(ops, in, block(ops, in), m)));
      });
}

void define_matrix(jlcxx::Module &mod) {
  using namespace xdiag;
  register_matrix<Spinhalf>(mod);
  register_matrix<tJ>(mod);
  register_matrix<Electron>(mod);
  register_matrix<Boson>(mod);
  register_matrix<Fermion>(mod);
}

} // namespace xdiag::julia
