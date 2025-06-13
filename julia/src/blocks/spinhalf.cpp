// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

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
      .method("nsites",
              [](Spinhalf const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
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

  mod.method("construct_Spinhalf",
             []() { JULIA_XDIAG_CALL_RETURN(Spinhalf()); });
  mod.method("construct_Spinhalf", [](int64_t nsites, std::string backend) {
    JULIA_XDIAG_CALL_RETURN(Spinhalf(nsites, backend));
  });
  mod.method("construct_Spinhalf",
             [](int64_t nsites, int64_t nup, std::string backend) {
               JULIA_XDIAG_CALL_RETURN(Spinhalf(nsites, nup, backend));
             });
  mod.method("construct_Spinhalf",
             [](int64_t nsites, Representation irrep, std::string backend) {
               JULIA_XDIAG_CALL_RETURN(Spinhalf(nsites, irrep, backend));
             });
  mod.method("construct_Spinhalf",
             [](int64_t nsites, int64_t nup, Representation irrep,
                std::string backend) {
               JULIA_XDIAG_CALL_RETURN(Spinhalf(nsites, nup, irrep, backend));
             });
}

} // namespace xdiag::julia
