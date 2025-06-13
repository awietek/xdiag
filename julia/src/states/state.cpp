// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "state.hpp"

namespace xdiag::julia {
void define_state(jlcxx::Module &mod) {

  mod.add_type<State>("cxx_State")
      .method("nsites",
              [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.nsites()) })
      .method("isreal",
              [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.isreal()) })
      .method("real", [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.real()) })
      .method("imag", [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.imag()) })
      .method("make_complex",
              [](State &s) { JULIA_XDIAG_CALL_VOID(s.make_complex()) })
      .method("dim", [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.dim()) })
      .method("size", [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.size()) })
      .method("nrows",
              [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.nrows()) })
      .method("ncols",
              [](State const &s) { JULIA_XDIAG_CALL_RETURN(s.ncols()) })
      .method("col", [](State const &s, int64_t n,
                        bool copy) { JULIA_XDIAG_CALL_RETURN(s.col(n, copy)) })
      .method("vector",
              [](State const &s, int64_t n, bool copy) {
                JULIA_XDIAG_CALL_RETURN(s.vector(n, copy))
              })
      .method("matrix",
              [](State const &s, bool copy) {
                JULIA_XDIAG_CALL_RETURN(s.matrix(copy))
              })
      .method("vectorC",
              [](State const &s, int64_t n, bool copy) {
                JULIA_XDIAG_CALL_RETURN(s.vectorC(n, copy))
              })
      .method("matrixC", [](State const &s, bool copy) {
        JULIA_XDIAG_CALL_RETURN(s.matrixC(copy))
      });

  mod.method("to_string", [](State const &s) { return to_string(s); });
  mod.method("isapprox",
             [](State const &v, State const &w, double rtol, double atol) {
               JULIA_XDIAG_CALL_RETURN(isapprox(v, w, rtol, atol));
             });

  mod.method("construct_State", []() { JULIA_XDIAG_CALL_RETURN(State()); });
  mod.method("construct_State",
             [](Spinhalf const &block, bool real, int64_t ncols) {
               JULIA_XDIAG_CALL_RETURN(State(block, real, ncols));
             });
  mod.method("construct_State", [](Spinhalf const &block, double const *ptr,
                                   int64_t ncols, int64_t stride) {
    JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols, stride));
  });
  mod.method("construct_State",
             [](Spinhalf const &block, complex const *ptr, int64_t ncols) {
               JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols));
             });
  mod.method("construct_State", [](tJ const &block, bool real, int64_t ncols) {
    JULIA_XDIAG_CALL_RETURN(State(block, real, ncols));
  });
  mod.method("construct_State", [](tJ const &block, double const *ptr,
                                   int64_t ncols, int64_t stride) {
    JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols, stride));
  });
  mod.method("construct_State",
             [](tJ const &block, complex const *ptr, int64_t ncols) {
               JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols));
             });
  mod.method("construct_State",
             [](Electron const &block, bool real, int64_t ncols) {
               JULIA_XDIAG_CALL_RETURN(State(block, real, ncols));
             });
  mod.method("construct_State", [](Electron const &block, double const *ptr,
                                   int64_t ncols, int64_t stride) {
    JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols, stride));
  });
  mod.method("construct_State",
             [](Electron const &block, complex const *ptr, int64_t ncols) {
               JULIA_XDIAG_CALL_RETURN(State(block, ptr, ncols));
             });
}

} // namespace xdiag::julia
