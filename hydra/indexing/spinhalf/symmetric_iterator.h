#pragma once

#include <hydra/common.h>

namespace hydra::indexing::spinhalf {

template <typename bit_t> class SymmetricIterator {
public:
  SymmetricIterator() = default;
  SymmetricIterator(std::vector<bit_t> const &reps, idx_t idx);

  inline bool operator==(SymmetricIterator<bit_t> const &rhs) const {
    return (idx_ == rhs.idx_) && (data_ == rhs.data_);
  }
  inline bool operator!=(SymmetricIterator<bit_t> const &rhs) const {
    return !operator==(rhs);
  }
  inline SymmetricIterator &operator++() {
    ++data_;
    ++idx_;
    return *this;
  }
    
  inline std::pair<bit_t, idx_t> operator*() const { return {*data_, idx_}; }

private:
  bit_t const *data_;
  idx_t idx_;
};

} // namespace hydra::indexing::spinhalf
