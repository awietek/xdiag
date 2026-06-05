// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <xdiag/extern/flat_hash_map.hpp>
#include <xdiag/extern/gsl/span>

#include <xdiag/basis/basis.hpp>
#include <xdiag/utils/likely.hpp>
#include <xdiag/utils/type_name.hpp>

#include <xdiag/symmetries/action/site_permutation_sublattice.hpp>
#include <xdiag/symmetries/permutation_group.hpp>
#include <xdiag/symmetries/representation.hpp>

namespace xdiag::basis {

constexpr int64_t maximum_prefix_bits = 32;

template <typename bit_tt, int n_sublat>
class BasisSublattice : public BasisType<BasisSublattice<bit_tt, n_sublat>> {
public:
  using bit_t = bit_tt;
  using iterator_t = typename std::vector<bit_t>::const_iterator;
  static constexpr std::string_view type_name =
      utils::get_type_name<BasisSublattice<bit_t, n_sublat>>();

  BasisSublattice() = default;
  BasisSublattice(PermutationGroup const &group, Vector const &characters);
  BasisSublattice(int64_t nup, PermutationGroup const &group,
                  Vector const &characters);

  iterator_t begin() const;
  iterator_t end() const;
  int64_t nsites() const;
  int64_t size() const override;
  int64_t d() const; // Local Hilbert space dimension per site
  int64_t dim() const;
  Vector const& characters() const;
  symmetries::SitePermutationSublattice<bit_t, n_sublat> const& action() const;

  int64_t index(bit_t state) const;
  inline bit_t operator[](int64_t idx) const { return reps_[idx]; }
  inline double norm(int64_t idx) const { return norms_[idx]; }
  inline double inv_norm(int64_t idx) const { return 1.0 / norms_[idx]; }

  ProductState
  product_state(int64_t idx,
                std::vector<std::string> const &dict) const override;

  bool operator==(BasisSublattice<bit_t, n_sublat> const &rhs) const;
  bool operator!=(BasisSublattice<bit_t, n_sublat> const &rhs) const;

private:
  PermutationGroup group_;
  symmetries::SitePermutationSublattice<bit_t, n_sublat> action_;
  Vector characters_;

  int64_t nsites_;
  int64_t nup_;
  int64_t n_postfix_bits_;

  std::vector<bit_t> reps_;
  std::vector<double> norms_;
  ska::flat_hash_map<bit_t, gsl::span<bit_t const>> rep_search_range_;

  int64_t index_of_representative(bit_t rep) const;

public:
  inline bit_t representative(bit_t state) const {
    return action_.representative(state);
  }

  // Returns {raw_rep_idx, sym, norm_out}; raw_rep_idx == 0 means zero-norm.
  // Actual index is raw_rep_idx - 1.
  inline std::tuple<int64_t, int64_t, double>
  representative_data(bit_t raw_state) const {
    auto [rep, sym] = action_.representative_sym(raw_state);
    int64_t idx = index_of_representative(rep);
    if (XDIAG_LIKELY(idx >= 0)) {
      return {idx + 1, sym, norms_[idx]};
    } else {
      return {0, 0, 0.0};
    }
  }

  std::pair<int64_t, int64_t> index_sym(bit_t raw_state) const;
  std::pair<int64_t, gsl::span<int64_t const>>
  index_syms(bit_t raw_state) const;
};

template <int n_sublat>
using BasisSublattice32 = BasisSublattice<uint32_t, n_sublat>;
template <int n_sublat>
using BasisSublattice64 = BasisSublattice<uint64_t, n_sublat>;
} // namespace xdiag::basis
