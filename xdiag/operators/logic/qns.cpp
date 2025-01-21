#include "qns.hpp"

#include <set>

#include <xdiag/operators/logic/isapprox.hpp>
#include <xdiag/operators/logic/permute.hpp>
#include <xdiag/operators/logic/valid.hpp>
#include <xdiag/utils/scalar.hpp>

namespace xdiag {

Spinhalf block(OpSum const &ops, Spinhalf const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    return (nup == nupr) ? block : Spinhalf(n_sites, *nupr, backend);
  } else if (!nup && irrep) {
    auto irrepr = representation(ops, block);
    return isapprox(irrep, irrepr) ? block : Spinhalf(n_sites, *irrep, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && isapprox(irrep, irrepr))
               ? block
               : Spinhalf(n_sites, *nupr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

tJ block(OpSum const &ops, tJ const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : tJ(n_sites, *nupr, *ndnr, backend);
  }
  // // Not yet implemented
  // else if (!nup && irrep) {
  //   auto irrepr = representation(ops, block);
  //   return isapprox(irrep, irrepr) ? block : tJ(n_sites, *irrep, backend);
  // }
  else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && (ndn == ndnr) && isapprox(irrep, irrepr))
               ? block
               : tJ(n_sites, *nupr, *ndnr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Electron block(OpSum const &ops, Electron const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  auto irrep = block.irrep();
  if (!nup && !irrep) {
    return block;
  } else if (nup && !irrep) {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : Electron(n_sites, *nupr, *ndnr, backend);
  } else if (!nup && irrep) {
    auto irrepr = representation(ops, block);
    return isapprox(irrep, irrepr) ? block : Electron(n_sites, *irrep, backend);
  } else { //(nup && irrep)
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    auto irrepr = representation(ops, block);
    return ((nup == nupr) && (ndn == ndnr) && isapprox(irrep, irrepr))
               ? block
               : Electron(n_sites, *nupr, *ndnr, *irrepr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

#ifdef XDIAG_USE_MPI
SpinhalfDistributed block(OpSum const &ops,
                          SpinhalfDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  if (!nup) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    return (nup == nupr) ? block : SpinhalfDistributed(n_sites, *nupr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
tJDistributed block(OpSum const &ops, tJDistributed const &block) try {
  int64_t n_sites = block.n_sites();
  std::string backend = block.backend();
  auto nup = block.nup();
  auto ndn = block.ndn();
  if (!nup) {
    return block;
  } else {
    auto nupr = nup(ops, block);
    auto ndnr = ndn(ops, block);
    return ((nup == nupr) && (ndn == ndnr))
               ? block
               : tJDistributed(n_sites, *nupr, *ndnr, backend);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
#endif

template <typename block_t>
Representation representation(OpSum const &ops, block_t const &block) try {
  auto irrep_block = block.irrep();
  if (!irrep_block) {
    XDIAG_THROW("Block has no irreducible representation defined. Hence, the "
                "irreducible representation of an OpSum times the block cannot "
                "be determined");
  }

  auto irrep_ops = representation(ops, irrep_block.group());
  if (irrep_ops) {
    return (*irrep_ops) * (*irrep_block);
  } else {
    XDIAG_THROW("OpSum does not transform according to an irredicuble "
                "representation of the symmetry group of the block")
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
#endif

template Representation representation(OpSum const &, Spinhalf const &);
template Representation representation(OpSum const &, tJ const &);
template Representation representation(OpSum const &, Electron const &);

Representation representation(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return representation(ops, b) },
                    block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

Representation representation(OpSum const &ops, State const &v) try {
  return representation(ops, state.block());
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
int64_t nup(OpSum const &ops, block_t const &block) try {
  auto nup_block = block.nup();
  if (!nup_block) {
    XDIAG_THROW("Block does not have a fixed number of nup particles");
  }

  auto nup_ops = nup(ops, irrep_block.group());
  if (nup_ops) {
    return (*nup_ops) + (*nup_block);
  } else {
    XDIAG_THROW("OpSum does not conserve the number of nup particles");
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
#endif

int64_t nup(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return nup(ops, b) }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t nup(OpSum const &ops, State const &v) try {
  return nup(ops, state.block());
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

template <typename block_t>
int64_t ndn(OpSum const &ops, block_t const &block) try {
  auto ndn_block = block.ndn();
  if (!ndn_block) {
    XDIAG_THROW("Block does not have a fixed number of ndn particles");
  }

  auto ndn_ops = ndn(ops, irrep_block.group());
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
#endif

int64_t ndn(OpSum const &ops, Block const &block) try {
  return std::visit([&](auto const &b) { return ndn(ops, b) }, block);
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}

int64_t ndn(OpSum const &ops, State const &v) try {
  return ndn(ops, state.block());
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
  check_valid(ops);
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> nups;
    for (auto [cpl, op] : ops) {
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
  check_valid(ops);
  if (ops.size() == 0) {
    return 0;
  } else {
    // Compute nup for every Op
    std::vector<std::optional<int64_t>> ndns;
    for (auto [cpl, op] : ops) {
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
  std::vector<complex> characters;
  std::vector<double> characters_real;
  bool real = true;
  for (auto const &perm : group) {
    OpSum opsp = permute(ops, perm);
    std::optional<Scalar> factor = isapprox_multiple(opsp, ops);
    if (factor) {
      characters.push_back((*factor).as<complex>());
      characters_real.push_back((*factor).real());
      if (!(*factor).isreal()) {
        real = false;
      }
    } else {
      return std::nullopt;
    }
  }
  if (real) {
    return Representation(group, characters_real);
  } else {
    return Representation(group, characters);
  }
} catch (Error const &e) {
  XDIAG_RETHROW(e);
}
} // namespace xdiag
