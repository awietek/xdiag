#include <julia/src/xdiagjl.hpp>

#include <julia/src/utils/armadillo.hpp>

#include <julia/src/operators/coupling.hpp>
#include <julia/src/operators/op.hpp>
#include <julia/src/operators/opsum.hpp>

#include <julia/src/blocks/electron.hpp>
#include <julia/src/blocks/spinhalf.hpp>
#include <julia/src/blocks/tj.hpp>

#include <julia/src/utils/utils.hpp>

#include <julia/src/states/state.hpp>
#include <julia/src/states/product_state.hpp>

// #include <julia/src/states/random_state.hpp>
// #include <julia/src/states/create_state.hpp>
// #include <julia/src/states/fill.hpp>

#include <julia/src/symmetries/permutation.hpp>
#include <julia/src/symmetries/permutation_group.hpp>
#include <julia/src/symmetries/representation.hpp>

#include <julia/src/algebra/matrix.hpp>

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

  // ProductState
  julia::define_product_state(mod);
  
  // Blocks
  julia::define_spinhalf(mod);
  julia::define_tj(mod);
  julia::define_electron(mod);

  // States
  julia::define_state(mod);
  
  // algebra
  julia::define_matrix(mod);

  // algorithms
  julia::define_eig0(mod);
  julia::define_eigval0(mod);

  // Utils
  julia::define_say_hello(mod);
  julia::define_set_verbosity(mod);

  // mod.add_type<State>("State")
  //     .constructor<>()
  //     .constructor<Spinhalf const &, bool, int64_t>()
  //     .constructor<Spinhalf const &, double *, int64_t>()
  //     .constructor<Spinhalf const &, complex *, int64_t>()
  //     .constructor<tJ const &, bool, int64_t>()
  //     .constructor<tJ const &, double *, int64_t>()
  //     .constructor<tJ const &, complex *, int64_t>()
  //     .constructor<Electron const &, bool, int64_t>()
  //     .constructor<Electron const &, double *, int64_t>()
  //     .constructor<Electron const &, complex *, int64_t>()
  //     .method("n_sites", &State::n_sites)
  //     .method("n_rows", &State::n_rows)
  //     .method("n_cols", &State::n_cols)
  //     .method("isreal", &State::isreal)
  //     .method("iscomplex", &State::iscomplex)
  //     .method("real", &State::real)
  //     .method("imag", &State::imag)
  //     .method("make_complex!", &State::make_complex);

  // mod.add_type<ProductState>("ProductStateCxx")
  //     .constructor<>()
  //     .constructor<std::vector<std::string> const &>()
  //     .method("n_sites", &ProductState::n_sites);

  // mod.method("zeros_like",
  //            [](State const &v) { JULIA_XDIAG_CALL_RETURN(zeros_like(v)); });

  // mod.method("memptr",
  //            [](State &state) { JULIA_XDIAG_CALL_RETURN(state.memptr()); });
  // mod.method("colptr", [](State &state, int64_t col) {
  //   JULIA_XDIAG_CALL_RETURN(state.colptr(col));
  // });

  // mod.method("memptrC",
  //            [](State &state) { JULIA_XDIAG_CALL_RETURN(state.memptrC()); });
  // mod.method("colptrC", [](State &state, int64_t col) {
  //   JULIA_XDIAG_CALL_RETURN(state.colptrC(col));
  // });

  // // methods to apply bonds
  // mod.method("apply_cxx", [](BondList const &bonds, State const &v, State &w)
  // {
  //   JULIA_XDIAG_CALL_VOID(apply(bonds, v, w));
  // });

  // // algebra operations
  // mod.method("norm_cxx",
  //            [](State const &v) { JULIA_XDIAG_CALL_RETURN(norm(v)); });

  // mod.method("dot_cxx", [](State const &v, State const &w) {
  //   JULIA_XDIAG_CALL_RETURN(dot(v, w));
  // });

  // mod.method("dotC_cxx", [](State const &v, State const &w) {
  //   JULIA_XDIAG_CALL_RETURN(dotC(v, w));
  // });

  // // Expectation values
  // mod.method("inner_cxx", [](Bond const &bond, State const &v) {
  //   JULIA_XDIAG_CALL_RETURN(inner(bond, v));
  // });

  // mod.method("innerC_cxx", [](Bond const &bond, State const &v) {
  //   JULIA_XDIAG_CALL_RETURN(innerC(bond, v));
  // });

  // mod.method("inner_cxx", [](State const &v, Bond const &bond, State const
  // &w) {
  //   JULIA_XDIAG_CALL_RETURN(inner(v, bond, w));
  // });

  // mod.method("innerC_cxx",
  //            [](State const &v, Bond const &bond, State const &w) {
  //              JULIA_XDIAG_CALL_RETURN(innerC(v, bond, w));
  //            });

  // mod.method("inner_cxx", [](BondList const &bonds, State const &v) {
  //   JULIA_XDIAG_CALL_RETURN(inner(bonds, v));
  // });

  // mod.method("innerC_cxx", [](BondList const &bonds, State const &v) {
  //   JULIA_XDIAG_CALL_RETURN(innerC(bonds, v));
  // });

  // mod.method("inner_cxx",
  //            [](State const &v, BondList const &bonds, State const &w) {
  //              JULIA_XDIAG_CALL_RETURN(inner(v, bonds, w));
  //            });

  // mod.method("innerC_cxx",
  //            [](State const &v, BondList const &bonds, State const &w) {
  //              JULIA_XDIAG_CALL_RETURN(innerC(v, bonds, w));
  //            });

}
