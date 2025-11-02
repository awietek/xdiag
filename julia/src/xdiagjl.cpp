// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#ifndef _OPENMP
#error "XDiag Julia wrapper needs to be compiled with OpenMP support"
#endif

#include <julia/src/xdiagjl.hpp>

#include <julia/src/utils/armadillo.hpp>

#include <julia/src/operators/block.hpp>
#include <julia/src/operators/hc.hpp>
#include <julia/src/operators/op.hpp>
#include <julia/src/operators/opsum.hpp>
#include <julia/src/operators/symmetrize.hpp>

#include <julia/src/blocks/electron.hpp>
#include <julia/src/blocks/spinhalf.hpp>
#include <julia/src/blocks/tj.hpp>

#include <julia/src/io/file_toml.hpp>
#include <julia/src/io/read.hpp>

#include <julia/src/utils/utils.hpp>

#include <julia/src/states/create_state.hpp>
#include <julia/src/states/fill.hpp>
#include <julia/src/states/gpwf.hpp>
#include <julia/src/states/product_state.hpp>
#include <julia/src/states/random_state.hpp>
#include <julia/src/states/state.hpp>

#include <julia/src/symmetries/permutation.hpp>
#include <julia/src/symmetries/permutation_group.hpp>
#include <julia/src/symmetries/representation.hpp>

#include <julia/src/algebra/algebra.hpp>
#include <julia/src/algebra/apply.hpp>
#include <julia/src/algebra/matrix.hpp>
#include <julia/src/algebra/sparse/coo_matrix.hpp>
#include <julia/src/algebra/sparse/csr_matrix.hpp>
#include <julia/src/algebra/sparse/csc_matrix.hpp>
#include <julia/src/algebra/sparse/apply.hpp>

#include <julia/src/algorithms/lanczos/eigs_lanczos.hpp>
#include <julia/src/algorithms/lanczos/eigvals_lanczos.hpp>
#include <julia/src/algorithms/sparse_diag.hpp>
#include <julia/src/algorithms/time_evolution/evolve_lanczos.hpp>
#include <julia/src/algorithms/time_evolution/imaginary_time_evolve.hpp>
#include <julia/src/algorithms/time_evolution/time_evolve.hpp>
#include <julia/src/algorithms/time_evolution/time_evolve_expokit.hpp>

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace xdiag;

  //
  // Comment: order of "defines" matters, otherwise errors like:
  // No appropriate factory for type ...
  //
  
  // Armadillo
  julia::define_vectors(mod);
  julia::define_matrices(mod);

  // Operators
  julia::define_op(mod);
  julia::define_opsum(mod);
  julia::define_hc(mod);

  // Symmetries
  julia::define_permutation(mod);
  julia::define_permutation_group(mod);
  julia::define_representation(mod);
  julia::define_symmetrize(mod);

  // ProductState
  julia::define_product_state(mod);

  // Blocks
  julia::define_spinhalf(mod);
  julia::define_tj(mod);
  julia::define_electron(mod);
  julia::define_block(mod);

  // States
  julia::define_state(mod);
  julia::define_random_state(mod);
  julia::define_gpwf(mod);
  julia::define_fill(mod);
  julia::define_create_state(mod);

  // Algebra
  julia::define_matrix(mod);
  julia::define_apply(mod);
  julia::define_algebra(mod);
  julia::define_coo_matrix(mod);
  julia::define_csr_matrix(mod);
  julia::define_csc_matrix(mod);
  julia::define_sparse_apply(mod);

  // Algorithms
  julia::define_eig0(mod);
  julia::define_eigs_lanczos(mod);
  julia::define_eigvals_lanczos(mod);
  julia::define_time_evolve(mod);
  julia::define_imaginary_time_evolve(mod);
  julia::define_time_evolve_expokit(mod);
  julia::define_evolve_lanczos(mod);

  // IO
  julia::define_file_toml(mod);
  julia::define_read(mod);

  // Utils
  julia::define_say_hello(mod);
  julia::define_print_version(mod);
  julia::define_set_verbosity(mod);
}
