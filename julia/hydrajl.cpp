#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"
#include "jlcxx/array.hpp"
#include "jlcxx/tuple.hpp"
#include "jlcxx/const_array.hpp"

#include <hydra/all.h>

// namespace hydra::julia {
// auto matrix(std::vector<std::string> const &bond_types,
//             std::vector<std::string> const &bond_couplings,
//             std::vector<int> const &bond_n_sites,
//             std::vector<int> const &bond_sites,
//             std::vector<std::string> const &coupling_name,
//             std::vector<std::string> const &coupling_value) {
//   static double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
//   return 1.2; // jlcxx::make_julia_array(d, 3, 2);
// }
// } // namespace hydra::julia


JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
  using namespace hydra;

  mod.add_type<Bond>("BondCxx")
      .constructor<>()
      .constructor<std::string, std::string, std::vector<int64_t> const &>();

  mod.add_type<BondList>("BondListCxx")
    .constructor<>()
    .constructor<std::vector<Bond> const&>()
    .method("set_coupling", &BondList::set_coupling);
  
  mod.method("matrix", []() {
    static double d[2][3] = {{1., 2., 3}, {4., 5., 6.}};
    return 1.2 ;//jlcxx::make_julia_array(d, 3, 2);
  });
  
  mod.add_type<Spinhalf>("Spinhalf")
    .constructor<int>()
    .constructor<int, int>()
    .method("n_sites", &Spinhalf::n_sites)
    .method("size", &Spinhalf::size);

}
