// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

namespace xdiag::julia {
void define_product_state(jlcxx::Module &mod) {

  using iterator_t = std::vector<std::string>::const_iterator;

  mod.add_type<iterator_t>("cxx_vector_string_iterator")
      .method("_incr", [](iterator_t &it) { ++it; })
      .method("_deref", [](iterator_t const &it) { return *it; })
      .method("==", [](iterator_t const &it1, iterator_t const &it2) {
        return it1 == it2;
      });

  mod.add_type<ProductState>("cxx_ProductState")
      .constructor<>()
      .constructor<int64_t>()
      .constructor<std::vector<std::string> const &>()
      .method(
          "nsites",
          [](ProductState const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("size",
              [](ProductState const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("getindex", [](ProductState const &s,
                             int64_t idx) { JULIA_XDIAG_CALL_RETURN(s[idx-1]) })
      .method("setindex!",
              [](ProductState &s, std::string local_state, int64_t idx) {
                JULIA_XDIAG_CALL_VOID(s[idx-1] = local_state)
              })
      .method("push!",
              [](ProductState &s, std::string local_state) {
                JULIA_XDIAG_CALL_VOID(s.push_back(local_state))
              })
      .method("_begin",
              [](ProductState const &s) { JULIA_XDIAG_CALL_RETURN(s.begin()) })
      .method("_end",
              [](ProductState const &s) { JULIA_XDIAG_CALL_RETURN(s.end()) });

  mod.method("to_string", [](ProductState const &s) { return to_string(s); });
}

} // namespace xdiag::julia
