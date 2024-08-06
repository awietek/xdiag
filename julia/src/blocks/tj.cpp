#include "tj.hpp"

namespace xdiag::julia {

void define_tj(jlcxx::Module &mod) {
  using iterator_t = tJIterator;
  mod.add_type<iterator_t>("cxx_tJIterator")
      .method("_incr", [](iterator_t &it) { ++it; })
      .method("_deref", [](iterator_t const &it) { return *it; })
      .method("!=", [](iterator_t const &it1, iterator_t const &it2) {
        return it1 != it2;
      });

  mod.add_type<tJ>("cxx_tJ")
      .constructor<int64_t, int64_t, int64_t>()
      .constructor<int64_t, int64_t, int64_t, PermutationGroup,
                   Representation>()
      .method("n_sites",
              [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.n_sites()) })
      .method("n_up", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.n_up()) })
      .method("n_dn", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.n_dn()) })
      .method(
          "permutation_group",
          [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.permutation_group()) })
      .method("irrep", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.irrep()) })
      .method("isreal", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) })
      .method("dim", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.dim()) })
      .method("size", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("_begin", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.begin()) })
      .method("_end", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.end()) })
      .method("index", [](tJ const &s, ProductState const &p) {
        JULIA_XDIAG_CALL_RETURN(s.index(p))
      });
  mod.method("to_string", [](tJ const &r) { return to_string(r); });
}

} // namespace xdiag::julia
