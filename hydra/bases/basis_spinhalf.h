// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef HYDRA_BASES_BASIS_SPINHALF_H_
#define HYDRA_BASES_BASIS_SPINHALF_H_

#include <string>

#include <hydra/common.h>
#include <hydra/utils/combinatorics.h>

#include <hydra/qns/qn_spinhalf.h>
#include <hydra/states/state_spinhalf.h>


namespace hydra {

template <class bit_t = std_bit_t> class BasisSpinHalfIterator {
public:
  using state_t = state_spinhalf<bit_t>;
  
  BasisSpinHalfIterator() = default;
  BasisSpinHalfIterator(const state_t &state);

  inline bool operator==(const BasisSpinHalfIterator<bit_t> &rhs) const {
    return current_ == rhs.current_;
  }
  inline bool operator!=(const BasisSpinHalfIterator<bit_t> &rhs) const {
    return !operator==(rhs);
  }
  inline BasisSpinHalfIterator &operator++() {
    current_ = combinatorics::get_next_pattern(current_);
    return *this;
  }
  inline state_t operator*() const { return state_t({current_}); }

private:
  bit_t current_;
};

template <class bit_t = std_bit_t> class BasisSpinHalf {
public:
  using qn_t = qn_spinhalf;
  using state_t = state_spinhalf<bit_t>;
  using iterator_t = BasisSpinHalfIterator<bit_t>;

  BasisSpinHalf() = default;
  BasisSpinHalf(int const &n_sites, qn_t const &qn);
  BasisSpinHalf(int const &n_sites, int const &qn);

  int n_sites() const;
  qn_t qn() const;
  iterator_t begin() const;
  iterator_t end() const;
  uint64 size() const;
  uint64 rawsize() const;

private:
  int n_sites_;
  qn_t qn_;
  iterator_t begin_, end_;
};

} // namespace hydra

#endif
