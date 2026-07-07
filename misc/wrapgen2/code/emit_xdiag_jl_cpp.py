#!/usr/bin/env python
from itertools import product
from common import *

def emit_xdiag_jl_cpp():
    total_string = """// SPDX-FileCopyrightText: 2026 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#ifndef _OPENMP
#error "XDiag Julia wrapper needs to be compiled with OpenMP support"
#endif

#include <julia/src/xdiagjl.hpp>
#include <julia/src/types.hpp>
#include <julia/src/modules.hpp>
    
JLCXX_MODULE define_julia_module(jlcxx::Module &mod) {
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
      .constructor<int64_t *, arma::uword, bool, bool>();

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
      .constructor<int64_t *, arma::uword, arma::uword, bool, bool>();
    
julia::define_types(mod);
julia::define_modules(mod);
}
"""
    return total_string

print(emit_xdiag_jl_cpp())
