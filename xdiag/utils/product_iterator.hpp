// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <utility>

namespace xdiag::utils {
template <typename I1, typename I2> class product_iterator {
public:
  product_iterator(I1 i1, I1 end1, I2 i2, I2 begin2, I2 end2)
      : i1_(i1), end1_(end1), i2_(i2), begin2_(begin2), end2_(end2) {}

  auto operator*() const { return std::make_pair(*i1_, *i2_); }

  product_iterator &operator++() {
    ++i2_;
    if (i2_ == end2_) {
      i2_ = begin2_;
      ++i1_;
    }
    return *this;
  }

  bool operator==(product_iterator const &other) const {
    return i1_ == other.i1_ && (i1_ == end1_ || i2_ == other.i2_);
  }

  bool operator!=(product_iterator const &other) const {
    return !(*this == other);
  }

private:
  I1 i1_, end1_;
  I2 i2_, begin2_, end2_;
};
} // namespace xdiag::utils
