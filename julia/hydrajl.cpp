#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <hydra/all.h>

#define JULIA_HYDRA_CALL_VOID(CMD)                                             \
  try {                                                                        \
    CMD;                                                                       \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in hydra C++ core library"));     \
  }

#define JULIA_HYDRA_CALL_RETURN(CMD)                                           \
  try {                                                                        \
    return CMD;                                                                \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in hydra C++ core library"));     \
    return decltype(CMD)();                                                    \
  }

#define JULIA_HYDRA_CALL_ASSIGN(LVALUE, CMD)                                   \
  try {                                                                        \
    LVALUE = CMD;                                                              \
  } catch (std::exception const &e) {                                          \
    traceback(e);                                                              \
    throw(std::runtime_error("Error occurred in hydra C++ core library"));     \
    return decltype(CMD)();                                                    \
  }

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace hydra;

  mod.add_type<Bond>("BondCxx")
      .constructor<>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>()
      .constructor<std::string, complex, std::vector<int64_t> const &>();

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

  mod.add_type<State>("State")
      .constructor<>()
      .constructor<Spinhalf const &, bool, int64_t>()
      .constructor<Spinhalf const &, double *, int64_t>()
      .constructor<Spinhalf const &, complex *, int64_t>()
      .method("n_sites", &State::n_sites)
      .method("n_rows", &State::n_rows)
      .method("n_cols", &State::n_cols)
      .method("isreal", &State::isreal)
      .method("iscomplex", &State::iscomplex)
      .method("real", &State::real)
      .method("imag", &State::imag)
      .method("make_complex!", &State::make_complex);

  mod.method("memptr",
             [](State &state) { JULIA_HYDRA_CALL_RETURN(state.memptr()); });
  mod.method("colptr", [](State &state, int64_t col) {
    JULIA_HYDRA_CALL_RETURN(state.colptr(col));
  });

  mod.method("memptrC",
             [](State &state) { JULIA_HYDRA_CALL_RETURN(state.memptrC()); });
  mod.method("colptrC", [](State &state, int64_t col) {
    JULIA_HYDRA_CALL_RETURN(state.colptrC(col));
  });

  // methods to compute matrices
  mod.method("matrix_cxx",
             [](double *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_HYDRA_CALL_VOID(matrix(mat, bonds, block_in, block_out));
             });

  mod.method("matrixC_cxx",
             [](complex *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               JULIA_HYDRA_CALL_VOID(matrixC(reinterpret_cast<complex *>(mat),
                                             bonds, block_in, block_out));
             });

  // Computing the ground state energy
  mod.method("eigval0_cxx", [](BondList const &bonds, Spinhalf const &block,
                               double precision, int max_iterations,
                               bool force_complex, uint64_t seed) {
    JULIA_HYDRA_CALL_RETURN(
        eigval0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  mod.method("eig0_cxx", [](BondList const &bonds, Spinhalf const &block,
                            double precision, int max_iterations,
                            bool force_complex, uint64_t seed) {
    JULIA_HYDRA_CALL_RETURN(
        eig0(bonds, block, precision, max_iterations, force_complex, seed));
  });

  // Print functions
  mod.method("print_pretty", [](const char *id, Spinhalf const &block) {
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
