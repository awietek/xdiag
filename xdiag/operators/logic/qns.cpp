// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "qns.hpp"

#include <optional>
#include <set>

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/operators/logic/order.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

template <typename block_t>
Representation representation(OpSum const &ops, block_t const &block) try {
  auto irrep_block = block.irrep();
  if (!irrep_block) {
    XDIAG_THROW("Block has no irreducible representation defined. Hence, the "
                "irreducible representation of an OpSum times the block cannot "
                "be determined");
  } else {
    auto irrep_ops = representation(ops, irrep_block->group());
    if (irrep_ops) {
      return (*irrep_ops) * (*irrep_block);
    } else {
      XDIAG_THROW("OpSum does not transform according to an irredicuble "
                  "representation of the symmetry group of the block")
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
template <>
Representation representation(OpSum const &ops,
                              SpinhalfDistributed const &block) try {
  XDIAG_THROW("Block of type SpinhalfDistributed does not have irreducible "
              "representation defined");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
Representation representation(OpSum const &ops,
                              tJDistributed const &block) try {
  XDIAG_THROW("Block of type tJDistributed does not have irreducible "
              "representation defined");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <>
Representation representation(OpSum const &ops,
                              ElectronDistributed const &block) try {
  XDIAG_THROW("Block of type ElectronDistributed does not have irreducible "
              "representation defined");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif

template Representation representation(OpSum const &, Spinhalf const &);
template Representation representation(OpSum const &, tJ const &);
template Representation representation(OpSum const &, Electron const &);

Representation representation(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return representation(ops, b); },
                    block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Representation representation(OpSum const &ops, State const &v) try {
  return representation(ops, v.block());
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
int64_t nup(OpSum const &ops, block_t const &block) try {
  auto nup_block = block.nup();
  if (!nup_block) {
    XDIAG_THROW("Block does not have a fixed number of nup particles");
  } else {
    auto nup_ops = nup(ops);
    if (nup_ops) {
      return (*nup_ops) + (*nup_block);
    } else {
      XDIAG_THROW("OpSum does not conserve the number of nup particles");
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template int64_t nup(OpSum const &ops, Spinhalf const &block);
template int64_t nup(OpSum const &ops, tJ const &block);
template int64_t nup(OpSum const &ops, Electron const &block);
#ifdef XDIAG_USE_MPI
template int64_t nup(OpSum const &ops, SpinhalfDistributed const &block);
template int64_t nup(OpSum const &ops, tJDistributed const &block);
template int64_t nup(OpSum const &ops, ElectronDistributed const &block);
#endif

int64_t nup(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return nup(ops, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t nup(OpSum const &ops, State const &v) try {
  return nup(ops, v.block());
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
int64_t ndn(OpSum const &ops, block_t const &block) try {
  auto ndn_block = block.ndn();
  if (!ndn_block) {
    XDIAG_THROW("Block does not have a fixed number of ndn particles");
  }

  auto ndn_ops = ndn(ops);
  if (ndn_ops) {
    return (*ndn_ops) + (*ndn_block);
  } else {
    XDIAG_THROW("OpSum does not conserve the number of ndn particles");
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <> int64_t ndn(OpSum const &ops, Spinhalf const &block) try {
  XDIAG_THROW("Block of type Spinhalf does not have ndn defined");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template int64_t ndn(OpSum const &ops, tJ const &block);
template int64_t ndn(OpSum const &ops, Electron const &block);
#ifdef XDIAG_USE_MPI
template <>
int64_t ndn(OpSum const &ops, SpinhalfDistributed const &block) try {
  XDIAG_THROW("Block of type SpinhalfDistributed does not have ndn defined");
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
template int64_t ndn(OpSum const &ops, tJDistributed const &block);
template int64_t ndn(OpSum const &ops, ElectronDistributed const &block);
#endif

int64_t ndn(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return ndn(ops, b); }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t ndn(OpSum const &ops, State const &v) try {
  return ndn(ops, v.block());
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename T>
static std::optional<int64_t> nup_matrix(arma::Mat<T> const &mat,
                                         double precision = 1e-12) {
  std::set<int64_t> diffs;
  int64_t diff;
  for (uint64_t i = 0; i < mat.n_rows; ++i) {
    for (uint64_t j = 0; j < mat.n_cols; ++j) {
      if (std::abs(mat(i, j)) > precision) {
        diff =
            bits::popcnt(i) - bits::popcnt(j); // needs generalization for d!=2
        diffs.insert(diff);
      }
    }
  }
  if (diffs.size() == 0) {
    return 0;
  } else if (diffs.size() == 1) {
    return diff;
  } else {
    return std::nullopt;
  }
}

std::optional<int64_t> nup(Op const &op) try {
  check_valid(op);
  std::string type = op.type();
  if ((type == "S+") || (type == "Cdagup")) {
    return 1;
  } else if ((type == "S-") || (type == "Cup")) {
    return -1;
  } else if (type == "Matrix") {
    Matrix mat = op.matrix();
    if (isreal(mat)) {
      return nup_matrix(mat.as<arma::mat>());
    } else {
      return nup_matrix(mat.as<arma::cx_mat>());
    }
  } else {
    return 0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> nup(OpSum const &ops) try {
  OpSum opso = order(ops);
  check_valid(opso);
  if (opso.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> nups;
    for (auto [cpl, op] : opso) {
      nups.push_back(nup(op));
    }

    // Check if all elements are the same
    if (std::adjacent_find(nups.begin(), nups.end(), std::not_equal_to<>()) ==
        nups.end()) {
      return nups[0];
    } else {
      return std::nullopt;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> ndn(Op const &op) try {
  check_valid(op);
  std::string type = op.type();
  if ((type == "S-") || (type == "Cdagdn")) {
    return 1;
  } else if ((type == "S+") || (type == "Cdn")) {
    return -1;
  } else {
    return 0;
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<int64_t> ndn(OpSum const &ops) try {
  OpSum opso = order(ops);
  check_valid(opso);
  if (opso.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> ndns;
    for (auto [cpl, op] : opso) {
      ndns.push_back(ndn(op));
    }

    // Check if all elements are the same
    if (std::adjacent_find(ndns.begin(), ndns.end(), std::not_equal_to<>()) ==
        ndns.end()) {
      return ndns[0];
    } else {
      return std::nullopt;
    }
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

std::optional<Representation>
representation(OpSum const &ops, PermutationGroup const &group) try {
  OpSum opso = order(ops);
  check_valid(opso);
  std::vector<complex> characters;
  bool real = true;
  for (auto const &perm : group) {
    OpSum opsp = permute(opso, perm);
    std::optional<Scalar> factor = isapprox_multiple(opsp, opso);
    if (factor) {
      characters.push_back((*factor).as<complex>());
    } else {
      return std::nullopt;
    }
  }
  return Representation(group, Vector(characters));
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag
