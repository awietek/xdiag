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

#ifndef HYDRA_BASES_BASIS_TJ_H_
#define HYDRA_BASES_BASIS_TJ_H_

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/qns/qn_tj.h>
#include <hydra/states/state_tj.h>

#include <iostream>
namespace hydra {

template <class bit_t = std_bit_t> class BasisTJIterator;

// BasisTJ
template <class bit_t = std_bit_t> class BasisTJ {
public:
  using qn_t = qn_tj;
  using state_t = state_tj<bit_t>;
  using iterator_t = BasisTJIterator<bit_t>;

  BasisTJ(int const &n_sites, qn_t const &qn);

  int n_sites() const;
  qn_t qn() const;
  iterator_t begin() const;
  iterator_t end() const;
  size_t size() const;
  size_t rawsize() const;

private:
  number_t n_sites_;
  qn_t qn_;
  iterator_t begin_, end_;
};

// BasisTJIterator
template <class bit_t> class BasisTJIterator {
public:
  BasisTJIterator() = default;
  BasisTJIterator(BasisSpinHalf<bit_t> const &up,
                  BasisSpinHalf<bit_t> const &holes,
                  state_tj<bit_t> const &state);

  inline bool operator==(BasisTJIterator<bit_t> const &rhs) const {
    return ((up_iter_ == rhs.up_iter_) && (holes_iter_ == rhs.holes_iter_));
  }

  inline bool operator!=(BasisTJIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }

  inline BasisTJIterator &operator++() {
    ++holes_iter_;
    if (holes_iter_ == holes_end_) {
      holes_iter_ = holes_begin_;
      ++up_iter_;
    }
    return *this;
  }

  inline state_tj<bit_t> operator*() const {
    bit_t ups = (*up_iter_).spins;
    bit_t dns = combinatorics::up_hole_to_down(ups, (*holes_iter_).spins);
    return {ups, dns};
  }

private:
  int n_sites_;
  BasisSpinHalfIterator<bit_t> holes_begin_, holes_end_;
  BasisSpinHalfIterator<bit_t> holes_iter_, up_iter_;
};

} // namespace hydra

#endif
