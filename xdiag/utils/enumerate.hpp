// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <tuple>
#include <utility>

namespace xdiag::utils {

template <typename container_t> class enumerate_wrapper_iterator {
public:
  enumerate_wrapper_iterator(int64_t index,
                             decltype(std::begin(std::declval<container_t &>()))
                                 const &iter)
      : index_(index), iter_(iter) {}
  inline bool operator!=(const iterator &other) const {
    return iter_ != other.iter_;
  }
  inline void operator++() {
    ++index_;
    ++iter_;
  }
  inline auto operator*() const { return std::tie(index_, *iter_); }

private:
  int64_t index_;
  decltype(std::begin(std::declval<container_t &>())) iter_;
};

template <typename container_t> class enumerate_wrapper {
public:
  using iterator_t = enumerate_wrapper_iterator<container_t>;
  explicit enumerate_wrapper(container_t &c) : container_(c) {}
  inline iterator_t begin() { return {0, std::begin(container_)}; }
  inline iterator_t end() { return {0, std::end(container_)}; }

private:
  container_t &container_;
};

template <typename container_t>
inline enumerate_wrapper<container_t> enumerate(container_t &c) {
  return enumerate_wrapper<container_t>(c);
}

} // namespace xdiag::utils
