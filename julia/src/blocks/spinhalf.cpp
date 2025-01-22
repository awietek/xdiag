#include "spinhalf.hpp"

namespace xdiag::julia {

void define_spinhalf(jlcxx::Module &mod) {
  using iterator_t = SpinhalfIterator;
  mod.add_type<iterator_t>("cxx_SpinhalfIterator")
      .method("_incr", [](iterator_t &it) { ++it; })
      .method("_deref", [](iterator_t const &it) { return *it; })
      .method("!=", [](iterator_t const &it1, iterator_t const &it2) {
        return it1 != it2;
      });

  mod.add_type<Spinhalf>("cxx_Spinhalf")
      .constructor<int64_t>()
      .constructor<int64_t, int64_t>()
      .constructor<int64_t, PermutationGroup, Representation>()
      .constructor<int64_t, int64_t, PermutationGroup, Representation>()
      .method("nsites",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("nup",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.nup()) })
      .method("permutation_group",
              [](Spinhalf const &s) {
                JULIA_XDIAG_CALL_RETURN(s.permutation_group())
              })
      .method("irrep",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.irrep()) })
      .method("isreal",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) })
      .method("dim", [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.dim()) })
      .method("size",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("_begin",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.begin()) })
      .method("_end",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.end()) })
      .method("index", [](Spinhalf const &s, ProductState const &p) {
        JULIA_XDIAG_CALL_RETURN(s.index(p))
      });

  mod.method("to_string", [](Spinhalf const &r) { return to_string(r); });
}

} // namespace xdiag::julia
