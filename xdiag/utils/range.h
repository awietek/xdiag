#pragma once

namespace xdiag::utils {

template <typename size_type = int> class range {
public:
  constexpr range(size_type e) : range(0, e) {}
  constexpr range(size_type b, size_type e) : current_(b), end_(e) {}
  size_type const &operator*() const { return current_; }
  range &operator++() {
    ++current_;
    return *this;
  }
  bool operator!=(const range &other) const {
    return current_ != other.current_;
  }

  range begin() const { return range(current_, end_); }
  range end() const { return range(end_, end_); }
private:
  size_type current_;
  size_type end_;
};

} // namespace xdiag::utils
