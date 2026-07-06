// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#ifndef _OPENMP
#error "XDiag Julia wrapper needs to be compiled with OpenMP support"
#endif

#include <julia/src/xdiagjl.hpp>
#include <julia/src/types.hpp>
#include <julia/src/modules.hpp>
// #include <julia/src/utils/armadillo.hpp>

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace xdiag;

  // // Armadillo bridge types (cxx_arma_vec/mat) must be registered before any
  // // generated method that returns or takes them.
  // julia::define_vectors(mod);
  // julia::define_matrices(mod);

  // All wrapped xdiag types + their methods (two-pass inside), then the
  // pointer-fill specials.
  julia::define_types(mod);
  julia::define_modules(mod);
}
