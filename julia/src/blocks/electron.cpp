#include "electron.hpp"

namespace xdiag::julia {

void define_electron(jlcxx::Module &mod) {

  using iterator_t = ElectronIterator;
  mod.add_type<iterator_t>("cxx_ElectronIterator")
      .method("_incr", [](iterator_t &it) { ++it; })
      .method("_deref", [](iterator_t const &it) { return *it; })
      .method("!=", [](iterator_t const &it1, iterator_t const &it2) {
        return it1 != it2;
      });

  mod.add_type<Electron>("cxx_Electron")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t, int64_t>()
      .constructor<int64_t, PermutationGroup, Representation>()
      .constructor<int64_t, int64_t, int64_t, PermutationGroup,
                   Representation>()
      .method("nsites",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("nup",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.nup()) })
      .method("ndn",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.ndn()) })
      .method("permutation_group",
              [](Electron const &s) {
                JULIA_XDIAG_CALL_RETURN(s.permutation_group())
              })
      .method("irrep",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.irrep()) })
      .method("isreal",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) })
      .method("dim", [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.dim()) })
      .method("size",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("_begin",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.begin()) })
      .method("_end",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.end()) })
      .method("index", [](Electron const &s, ProductState const &p) {
        JULIA_XDIAG_CALL_RETURN(s.index(p))
      });

  mod.method("to_string", [](Electron const &r) { return to_string(r); });
}
} // namespace xdiag::julia
