#include <julia/src/xdiagjl.hpp>

#include <julia/src/utils/armadillo.hpp>

#include <julia/src/operators/coupling.hpp>
#include <julia/src/operators/op.hpp>
#include <julia/src/operators/opsum.hpp>
#include <julia/src/operators/symmetrize.hpp>

#include <julia/src/blocks/electron.hpp>
#include <julia/src/blocks/spinhalf.hpp>
#include <julia/src/blocks/tj.hpp>

#include <julia/src/utils/utils.hpp>

#include <julia/src/states/state.hpp>
#include <julia/src/states/product_state.hpp>
#include <julia/src/states/random_state.hpp>
#include <julia/src/states/fill.hpp>
#include <julia/src/states/gpwf.hpp>
#include <julia/src/states/create_state.hpp>

#include <julia/src/symmetries/permutation.hpp>
#include <julia/src/symmetries/permutation_group.hpp>
#include <julia/src/symmetries/representation.hpp>

#include <julia/src/algebra/matrix.hpp>
#include <julia/src/algebra/apply.hpp>
#include <julia/src/algebra/algebra.hpp>

#include <julia/src/algorithms/sparse_diag.hpp>

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace xdiag;

  // Armadillo
  julia::define_vectors(mod);
  julia::define_matrices(mod);

  // Operators
  julia::define_coupling(mod);
  julia::define_op(mod);
  julia::define_opsum(mod);

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

  // States
  julia::define_state(mod);
  julia::define_random_state(mod);
  julia::define_gpwf(mod);
  julia::define_fill(mod);
  julia::define_create_state(mod);
  
  // algebra
  julia::define_matrix(mod);
  julia::define_apply(mod);
  julia::define_algebra(mod);

  // algorithms
  julia::define_eig0(mod);
  julia::define_eigval0(mod);

  // Utils
  julia::define_say_hello(mod);
  julia::define_print_version(mod);
  julia::define_set_verbosity(mod);

}
