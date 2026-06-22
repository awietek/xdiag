// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once
#ifdef XDIAG_DISTRIBUTED

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include <xdiag/armadillo.hpp>

#include <xdiag/algebra/algebra.hpp>
#include <xdiag/algebra/normal_order.hpp>
#include <xdiag/mpi/buffer.hpp>
#include <xdiag/operators/coeff.hpp>
#include <xdiag/operators/op.hpp>
#include <xdiag/operators/opsum.hpp>
#include <xdiag/utils/error.hpp>

#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_exchange.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_hopping.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_number.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_number_number.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_raise_lower.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_szsz.hpp>
#include <xdiag/kernels/blocks/distributed/electron_distributed/terms/apply_u.hpp>

namespace xdiag::basis::electron_distributed {

// Up-species term: hops/raises in the up spin (Hopup / Cdagup / Cup).
inline bool is_up_term(std::string const &type) {
  return (type == "Hopup") || (type == "Cdagup") || (type == "Cup");
}

// Matrix-free, MPI-aware application of a (number conserving in both species,
// no symmetry) OpSum on a distributed electron basis (double occupancy
// allowed). Same scheme as the t-J block: dn-species + diagonal terms act in
// the native up/dn ordering; up-species terms act after a transpose to dn/up
// ordering. Ops are normal-ordered to single elementary operators first.
template <typename coeff_t, class basis_t>
void apply_terms(OpSum const &ops, basis_t const &basis_in,
                 arma::Col<coeff_t> const &vec_in, basis_t const &basis_out,
                 arma::Col<coeff_t> &vec_out) try {

  auto algebra = algebra::electron_implementation_algebra(basis_in.nsites());
  auto ops_compiled = normal_order(ops.plain(), algebra);

  std::vector<std::pair<Coeff, Op>> terms;
  bool has_up = false;
  for (auto const &[c, monomial] : ops_compiled) {
    if (monomial.size() != 1) {
      XDIAG_THROW("ElectronDistributed only supports single-operator terms "
                  "(monomials of length 1). Products of operators are not "
                  "implemented.");
    }
    Op op = monomial[0];
    terms.push_back({c, op});
    has_up |= is_up_term(op.type());
  }

  int64_t buffer_size =
      std::max({basis_in.size(), basis_in.size_transpose(), basis_out.size(),
                basis_out.size_transpose()});
  mpi::buffer.reserve<coeff_t>(buffer_size);

  // Operators in the native up/dn ordering (dn species + diagonal terms).
  for (auto const &[c, op] : terms) {
    std::string type = op.type();
    if (type == "SzSz") {
      apply_szsz<coeff_t>(c, op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if ((type == "Nup") || (type == "Ndn")) {
      apply_number<coeff_t>(c, op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "Nupdn") {
      apply_nupdn<coeff_t>(c, op, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "NupdnNupdn") {
      apply_nupdn_nupdn<coeff_t>(c, op, basis_in, vec_in.memptr(),
                                 vec_out.memptr());
    } else if (type == "NtotNtot") {
      apply_ntot_ntot<coeff_t>(c, op, basis_in, vec_in.memptr(),
                               vec_out.memptr());
    } else if (type == "NupNdn") {
      apply_nup_ndn<coeff_t>(c, op, basis_in, vec_in.memptr(),
                             vec_out.memptr());
    } else if (type == "NupNup") {
      apply_nup_nup<coeff_t>(c, op, basis_in, vec_in.memptr(),
                             vec_out.memptr());
    } else if (type == "NdnNdn") {
      apply_ndn_ndn<coeff_t>(c, op, basis_in, vec_in.memptr(),
                             vec_out.memptr());
    } else if (type == "Exchange") {
      apply_exchange<coeff_t>(c, op, basis_in, vec_in.memptr(),
                              vec_out.memptr());
    } else if (type == "Hopdn") {
      apply_hopping<coeff_t>(c, op, basis_in, vec_in.memptr(),
                             vec_out.memptr());
    } else if ((type == "Cdagdn") || (type == "Cdn")) {
      apply_raise_lower<coeff_t>(c, op, basis_in, vec_in.memptr(), basis_out,
                                 vec_out.memptr());
    } else if (type == "HubbardU") {
      apply_u<coeff_t>(c, basis_in, vec_in.memptr(), vec_out.memptr());
    } else if (type == "Id") {
      coeff_t cc = c.scalar().as<coeff_t>();
      for (int64_t i = 0; i < basis_in.size(); ++i) {
        vec_out[i] += cc * vec_in[i];
      }
    } else if (is_up_term(type)) {
      continue; // handled below in the transposed ordering
    } else {
      XDIAG_THROW(
          std::string("Unsupported Op type for ElectronDistributed block: ") +
          type);
    }
  }

  // Up-species operators: transpose to dn/up ordering, apply, transpose back.
  if (has_up) {
    basis_in.transpose(vec_in.memptr());
    coeff_t *vec_in_trans = mpi::buffer.send<coeff_t>();
    coeff_t *vec_out_trans = mpi::buffer.recv<coeff_t>();

    for (auto const &[c, op] : terms) {
      std::string type = op.type();
      if (type == "Hopup") {
        apply_hopping<coeff_t>(c, op, basis_in, vec_in_trans, vec_out_trans);
      } else if ((type == "Cdagup") || (type == "Cup")) {
        apply_raise_lower<coeff_t>(c, op, basis_in, vec_in_trans, basis_out,
                                   vec_out_trans);
      }
    }

    basis_out.transpose_r(vec_out_trans);
    coeff_t *send = mpi::buffer.send<coeff_t>();
    for (int64_t i = 0; i < basis_out.size(); ++i) {
      vec_out[i] += send[i];
    }
  }
}
XDIAG_CATCH

} // namespace xdiag::basis::electron_distributed
#endif
