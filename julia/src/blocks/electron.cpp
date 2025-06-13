// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

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
      .method("nsites",
              [](Electron const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
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

  mod.method("construct_Electron", []() {
    JULIA_XDIAG_CALL_RETURN(Electron());
  });
  mod.method("construct_Electron", [](int64_t nsites, std::string backend) {
    JULIA_XDIAG_CALL_RETURN(Electron(nsites, backend));
  });
  mod.method("construct_Electron",
             [](int64_t nsites, int64_t nup, int64_t ndn, std::string backend) {
               JULIA_XDIAG_CALL_RETURN(Electron(nsites, nup, ndn, backend));
             });
  mod.method("construct_Electron",
             [](int64_t nsites, Representation irrep, std::string backend) {
               JULIA_XDIAG_CALL_RETURN(Electron(nsites, irrep, backend));
             });
  mod.method("construct_Electron", [](int64_t nsites, int64_t nup, int64_t ndn,
                                      Representation irrep,
                                      std::string backend) {
    JULIA_XDIAG_CALL_RETURN(Electron(nsites, nup, ndn, irrep, backend));
  });
}
} // namespace xdiag::julia
