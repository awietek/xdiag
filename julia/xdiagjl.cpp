#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <xdiag/all.hpp>

#define JULIA_XDIAG_CALL_VOID(CMD)                                             \
  try {                                                                        \
    CMD;                                                                       \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in xdiag C++ core library"));     \
  }

#define JULIA_XDIAG_CALL_RETURN(CMD)                                           \
  try {                                                                        \
    return CMD;                                                                \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in xdiag C++ core library"));     \
    return decltype(CMD)();                                                    \
  }

#define JULIA_XDIAG_CALL_ASSIGN(LVALUE, CMD)                                   \
  try {                                                                        \
    LVALUE = CMD;                                                              \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in xdiag C++ core library"));     \
    return decltype(CMD)();                                                    \
  }

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace xdiag;

  mod.add_type<Bond>("BondCxx")
      .constructor<>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>()
      .constructor<std::string, complex, std::vector<int64_t> const &>()
      .constructor<std::string, double, std::vector<int64_t> const &>()
      .method("isreal", &Bond::isreal)
      .method("iscomplex", &Bond::iscomplex);

  mod.add_type<BondList>("BondListCxx")
      .constructor<>()
      .constructor<std::vector<Bond> const &>()
      .method("set_coupling", &BondList::set_coupling)
      .method("isreal", &BondList::isreal)
      .method("iscomplex", &BondList::iscomplex);

  mod.add_type<Spinhalf>("Spinhalf")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t>()
      .method("n_sites", &Spinhalf::n_sites)
      .method("size", &Spinhalf::size)
      .method("isreal", &Spinhalf::isreal)
      .method("iscomplex", &Spinhalf::iscomplex);

  mod.add_type<tJ>("tJ")
      .constructor<int64_t, int64_t, int64_t>()
      .method("n_sites", &tJ::n_sites)
      .method("n_up", &tJ::n_up)
      .method("n_dn", &tJ::n_dn)
      .method("size", &tJ::size)
      .method("isreal", &tJ::isreal)
      .method("iscomplex", &tJ::iscomplex);

  mod.add_type<Electron>("Electron")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t, int64_t>()
      .method("n_sites", &Electron::n_sites)
      .method("n_up", &Electron::n_up)
      .method("n_dn", &Electron::n_dn)
      .method("size", &Electron::size)
      .method("isreal", &Electron::isreal)
      .method("iscomplex", &Electron::iscomplex);

  mod.add_type<State>("State")
      .constructor<>()
      .constructor<Spinhalf const &, bool, int64_t>()
      .constructor<Spinhalf const &, double *, int64_t>()
      .constructor<Spinhalf const &, complex *, int64_t>()
      .constructor<tJ const &, bool, int64_t>()
      .constructor<tJ const &, double *, int64_t>()
      .constructor<tJ const &, complex *, int64_t>()
      .constructor<Electron const &, bool, int64_t>()
      .constructor<Electron const &, double *, int64_t>()
      .constructor<Electron const &, complex *, int64_t>()
      .method("n_sites", &State::n_sites)
      .method("n_rows", &State::n_rows)
      .method("n_cols", &State::n_cols)
      .method("isreal", &State::isreal)
      .method("iscomplex", &State::iscomplex)
      .method("real", &State::real)
      .method("imag", &State::imag)
      .method("make_complex!", &State::make_complex);

  mod.add_type<ProductState>("ProductStateCxx")
      .constructor<>()
      .constructor<std::vector<std::string> const &>()
      .method("n_sites", &ProductState::n_sites);

  mod.method("zeros_like",
             [](State const &v) { JULIA_XDIAG_CALL_RETURN(zeros_like(v)); });

  mod.method("memptr",
             [](State &state) { JULIA_XDIAG_CALL_RETURN(state.memptr()); });
  mod.method("colptr", [](State &state, int64_t col) {
    JULIA_XDIAG_CALL_RETURN(state.colptr(col));
  });

  mod.method("memptrC",
             [](State &state) { JULIA_XDIAG_CALL_RETURN(state.memptrC()); });
  mod.method("colptrC", [](State &state, int64_t col) {
    JULIA_XDIAG_CALL_RETURN(state.colptrC(col));
  });

  // methods to compute matrices
  mod.method("matrix_cxx",
             [](double *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
             });

  mod.method("matrixC_cxx",
             [](complex *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrixC(reinterpret_cast<complex *>(mat),
                                             bonds, block_in, block_out));
             });

  mod.method("matrix_cxx", [](double *mat, BondList const &bonds,
                              tJ const &block_in, tJ const &block_out) {
    JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
  });

  mod.method("matrixC_cxx", [](complex *mat, BondList const &bonds,
                               tJ const &block_in, tJ const &block_out) {
    JULIA_XDIAG_CALL_VOID(
        matrixC(reinterpret_cast<complex *>(mat), bonds, block_in, block_out));
  });

  mod.method("matrix_cxx",
             [](double *mat, BondList const &bonds, Electron const &block_in,
                Electron const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrix(mat, bonds, block_in, block_out));
             });

  mod.method("matrixC_cxx",
             [](complex *mat, BondList const &bonds, Electron const &block_in,
                Electron const &block_out) {
               JULIA_XDIAG_CALL_VOID(matrixC(reinterpret_cast<complex *>(mat),
                                             bonds, block_in, block_out));
             });

  // methods to apply bonds
  mod.method("apply_cxx", [](BondList const &bonds, State const &v, State &w) {
    JULIA_XDIAG_CALL_VOID(apply(bonds, v, w));
  });

  // algebra operations
  mod.method("norm_cxx",
             [](State const &v) { JULIA_XDIAG_CALL_RETURN(norm(v)); });

  mod.method("dot_cxx", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(dot(v, w));
  });

  mod.method("dotC_cxx", [](State const &v, State const &w) {
    JULIA_XDIAG_CALL_RETURN(dotC(v, w));
  });

  // Expectation values
  mod.method("inner_cxx", [](Bond const &bond, State const &v) {
    JULIA_XDIAG_CALL_RETURN(inner(bond, v));
  });

  mod.method("innerC_cxx", [](Bond const &bond, State const &v) {
    JULIA_XDIAG_CALL_RETURN(innerC(bond, v));
  });

  mod.method("inner_cxx", [](State const &v, Bond const &bond, State const &w) {
    JULIA_XDIAG_CALL_RETURN(inner(v, bond, w));
  });

  mod.method("innerC_cxx",
             [](State const &v, Bond const &bond, State const &w) {
               JULIA_XDIAG_CALL_RETURN(innerC(v, bond, w));
             });

  mod.method("inner_cxx", [](BondList const &bonds, State const &v) {
    JULIA_XDIAG_CALL_RETURN(inner(bonds, v));
  });

  mod.method("innerC_cxx", [](BondList const &bonds, State const &v) {
    JULIA_XDIAG_CALL_RETURN(innerC(bonds, v));
  });

  mod.method("inner_cxx",
             [](State const &v, BondList const &bonds, State const &w) {
               JULIA_XDIAG_CALL_RETURN(inner(v, bonds, w));
             });

  mod.method("innerC_cxx",
             [](State const &v, BondList const &bonds, State const &w) {
               JULIA_XDIAG_CALL_RETURN(innerC(v, bonds, w));
             });

  // methods for time evolution
  mod.method("exp_sym_v_cxx", [](BondList const &bonds, State state, double tau,
                                 bool normalize, double shift, double precision,
                                 int64_t max_iterations, double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(exp_sym_v(bonds, state, tau, normalize, shift,
                                      precision, max_iterations,
                                      deflation_tol));
  });

  mod.method("exp_sym_v_cxx", [](BondList const &bonds, State state,
                                 complex tau, bool normalize, double shift,
                                 double precision, int64_t max_iterations,
                                 double deflation_tol) {
    JULIA_XDIAG_CALL_RETURN(exp_sym_v(bonds, state, tau, normalize, shift,
                                      precision, max_iterations,
                                      deflation_tol));
  });

  // Computing the ground state energy
  mod.method("eigval0_cxx", [](BondList const &bonds, Spinhalf const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eigval0_cxx", [](BondList const &bonds, tJ const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eigval0_cxx", [](BondList const &bonds, Electron const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0_cxx", [](BondList const &bonds, Spinhalf const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0_cxx", [](BondList const &bonds, tJ const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0_cxx", [](BondList const &bonds, Electron const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_XDIAG_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  // Print functions
  mod.method("print_pretty", [](const char *id, Spinhalf const &block) {
    utils::print_pretty(id, block);
  });

  mod.method("print_pretty", [](const char *id, tJ const &block) {
    utils::print_pretty(id, block);
  });

  mod.method("print_pretty", [](const char *id, Electron const &block) {
    utils::print_pretty(id, block);
  });

  mod.method("print_pretty", [](const char *id, State const &state) {
    utils::print_pretty(id, state);
  });

  mod.method("print_pretty", [](const char *id, Bond const &bond) {
    utils::print_pretty(id, bond);
  });

  mod.method("print_pretty", [](const char *id, BondList const &bonds) {
    utils::print_pretty(id, bonds);
  });
}
