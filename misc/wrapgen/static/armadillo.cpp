// SPDX-License-Identifier: Apache-2.0
// STATIC hand-written bridge -- copied verbatim by generate.sh.

#include <julia/src/armadillo.hpp>

namespace xdiag::julia {

void define_armadillo(jlcxx::Module &mod) {
  using namespace xdiag;

  mod.add_type<arma::vec>("typ_arma_vec")
      .constructor<double *, arma::uword, bool, bool>()
      .method("memptr", [](arma::vec &m) { return m.memptr(); })
      .method("n_rows", [](arma::vec const &m) { return m.n_rows; });

  mod.add_type<arma::cx_vec>("typ_arma_cx_vec")
      .constructor<complex *, arma::uword, bool, bool>()
      .method("memptr", [](arma::cx_vec &m) { return m.memptr(); })
      .method("n_rows", [](arma::cx_vec const &m) { return m.n_rows; });

  mod.add_type<arma::Col<int64_t>>("typ_arma_vec_int64_t")
      .constructor<int64_t *, arma::uword, bool, bool>()
      .method("memptr", [](arma::Col<int64_t> &m) { return m.memptr(); })
      .method("n_rows", [](arma::Col<int64_t> const &m) { return m.n_rows; });

  mod.add_type<arma::mat>("typ_arma_mat")
      .constructor<double *, arma::uword, arma::uword, bool, bool>()
      .method("memptr", [](arma::mat &m) { return m.memptr(); })
      .method("n_rows", [](arma::mat const &m) { return m.n_rows; })
      .method("n_cols", [](arma::mat const &m) { return m.n_cols; })
      .method("n_elem", [](arma::mat const &m) { return m.n_elem; });

  mod.add_type<arma::cx_mat>("typ_arma_cx_mat")
      .constructor<complex *, arma::uword, arma::uword, bool, bool>()
      .method("memptr", [](arma::cx_mat &m) { return m.memptr(); })
      .method("n_rows", [](arma::cx_mat const &m) { return m.n_rows; })
      .method("n_cols", [](arma::cx_mat const &m) { return m.n_cols; })
      .method("n_elem", [](arma::cx_mat const &m) { return m.n_elem; });

  mod.add_type<arma::Mat<int64_t>>("typ_arma_mat_int64_t")
      .constructor<int64_t *, arma::uword, arma::uword, bool, bool>()
      .method("memptr", [](arma::Mat<int64_t> &m) { return m.memptr(); })
      .method("n_rows", [](arma::Mat<int64_t> const &m) { return m.n_rows; })
      .method("n_cols", [](arma::Mat<int64_t> const &m) { return m.n_cols; })
      .method("n_elem", [](arma::Mat<int64_t> const &m) { return m.n_elem; });
}

} // namespace xdiag::julia
