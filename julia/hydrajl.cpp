#include "jlcxx/array.hpp"
#include "jlcxx/const_array.hpp"
#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/tuple.hpp"

#include <hydra/all.h>

JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace hydra;

  mod.add_type<Bond>("BondCxx")
      .constructor<>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>();

  mod.add_type<BondList>("BondListCxx")
      .constructor<>()
      .constructor<std::vector<Bond> const &>()
      .method("set_coupling", &BondList::set_coupling);

  mod.add_type<Spinhalf>("Spinhalf")
      .constructor<int>()
      .constructor<int, int>()
      .method("n_sites", &Spinhalf::n_sites)
      .method("size", &Spinhalf::size);

  // methods to compute matrices
  mod.method("matrix_real_cxx",
             [](double *mat, BondList const &bonds, Spinhalf const &block_in,
                Spinhalf const &block_out) {
               matrix_real(mat, bonds, block_in, block_out);
             });

  mod.method("matrix_cplx_cxx", [](double *mat, BondList const &bonds,
                                   Spinhalf const &block_in,
                                   Spinhalf const &block_out) {
    matrix_cplx(reinterpret_cast<complex *>(mat), bonds, block_in, block_out);
  });

  // Computing the ground state energy
  mod.method("eig0_real_cxx",
             [](BondList const &bonds, Spinhalf const &block, double precision,
                int max_iterations, uint64_t seed) {
               return eig0_real(bonds, block, precision, max_iterations, seed);
             });
  mod.method("eig0_cplx_cxx",
             [](BondList const &bonds, Spinhalf const &block, double precision,
                int max_iterations, uint64_t seed) {
               return eig0_cplx(bonds, block, precision, max_iterations, seed);
             });

  // Computing the ground state
  mod.method("groundstate_real_cxx", [](BondList const &bonds,
                                        Spinhalf const &block, double precision,
                                        int max_iterations, uint64_t seed) {
    return groundstate_real(bonds, block, precision, max_iterations, seed);
  });
  mod.method("groundstate_cplx_cxx", [](BondList const &bonds,
                                        Spinhalf const &block, double precision,
                                        int max_iterations, uint64_t seed) {
    return groundstate_cplx(bonds, block, precision, max_iterations, seed);
  });

  // Print functions
  mod.method("print_pretty", [](const char *id, Spinhalf const &block) {
    utils::print_pretty(id, block);
  });
}
