// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

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
      .constructor<>()
      .constructor<int64_t, int64_t, int64_t, std::string>()
      .constructor<int64_t, int64_t, int64_t, Representation, std::string>()
      .method("nsites", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("isreal", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) })
      .method("dim", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.dim()) })
      .method("size", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("_begin", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.begin()) })
      .method("_end", [](tJ const &s) { JULIA_XDIAG_CALL_RETURN(s.end()) })
      .method("index", [](tJ const &s, ProductState const &p) {
        JULIA_XDIAG_CALL_RETURN(s.index(p))
      });
  mod.method("to_string", [](tJ const &r) { return to_string(r); });

  mod.method("construct_tJ", []() {
    JULIA_XDIAG_CALL_RETURN(tJ());
  });
  mod.method("construct_tJ",
             [](int64_t nsites, int64_t nup, int64_t ndn, std::string backend) {
               JULIA_XDIAG_CALL_RETURN(tJ(nsites, nup, ndn, backend));
             });
  mod.method("construct_tJ", [](int64_t nsites, int64_t nup, int64_t ndn,
                                      Representation irrep,
                                      std::string backend) {
    JULIA_XDIAG_CALL_RETURN(tJ(nsites, nup, ndn, irrep, backend));
  });
}

} // namespace xdiag::julia
