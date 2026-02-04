// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0
#include <vector>
#include <xdiag/bits/bitvector.hpp>
#include <xdiag/symmetries/group_action/group_action.hpp>

namespace xdiag::symmetries {

template <typename state_iterator_tt> class RepresentativeLookupTable {
public:
  using state_iterator_t = state_iterator_tt;
  using bit_t = typename state_iterator_t::bit_t;

  RepresentativeLookupTable() = default;
  RepresentativeLookupTable(state_iterator_t const &state_iterator,
                            GroupAction const &group_action,
                            std::vector<double> const &characters);
  RepresentativeLookupTable(state_iterator_t const &state_iterator,
                            GroupAction const &group_action,
                            std::vector<complex> const &characters);

  inline int64_t nrepresentatives() const { return representative_.size(); }
  inline int64_t representative_index(int64_t idx) const {
    return representative_index_[idx];
  }
  inline int64_t representative_symmetry(int64_t idx) const {
    return representative_symmetry_[idx];
  }
  inline int64_t representative_norm(int64_t idx) const {
    return norm_[representative_norm_index_[idx]];
  }

  bool operator==(RepresentativeLookupTable<state_iterator_t> const &rhs) const;
  bool operator!=(RepresentativeLookupTable<state_iterator_t> const &rhs) const;

private:
  bits::BitVector<bit_t> representative_;
  bits::BitVector<bit_t> representative_index_;
  bits::BitVector<bit_t> representative_symmetry_;
  bits::BitVector<bit_t> representative_norm_index_;
  std::vector<double> norm_;
};

} // namespace xdiag::symmetries
