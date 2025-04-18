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
      .constructor<>()
      .constructor<int64_t, std::string>()
      .constructor<int64_t, int64_t, std::string>()
      .constructor<int64_t, Representation, std::string>()
      .constructor<int64_t, int64_t, Representation, std::string>()
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
}

} // namespace xdiag::julia
