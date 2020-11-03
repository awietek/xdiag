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

#ifndef HYDRA_BASES_BASIS_ELECTRON_H_
#define HYDRA_BASES_BASIS_ELECTRON_H_

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/qns/qn_electron.h>
#include <hydra/states/state_electron.h>

namespace hydra {

template <class bit_t = std_bit_t> class BasisElectronIterator;

// BasisElectron
template <class bit_t = std_bit_t> class BasisElectron {
public:
  using qn_t = qn_electron;
  using state_t = state_electron<bit_t>;
  using iterator_t = BasisElectronIterator<bit_t>;

  BasisElectron(int const &n_sites, qn_t const &qn);

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

// BasisElectronIterator
template <class bit_t> class BasisElectronIterator {
public:
  BasisElectronIterator() = default;
  BasisElectronIterator(BasisSpinHalf<bit_t> const &up,
                        BasisSpinHalf<bit_t> const &down,
                        state_electron<bit_t> const &state);

  inline bool operator==(BasisElectronIterator<bit_t> const &rhs) const {
    return ((up_iter_ == rhs.up_iter_) && (down_iter_ == rhs.down_iter_));
  }

  inline bool operator!=(BasisElectronIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }

  inline BasisElectronIterator &operator++() {
    ++down_iter_;
    if (down_iter_ == down_end_) {
      down_iter_ = down_begin_;
      ++up_iter_;
    }
    return *this;
  }

  inline state_electron<bit_t> operator*() const {
    return {(*up_iter_).spins, (*down_iter_).spins};
  }

private:
  int n_sites_;
  BasisSpinHalfIterator<bit_t> down_begin_, down_end_;
  BasisSpinHalfIterator<bit_t> down_iter_, up_iter_;
};

} // namespace hydra

#endif
