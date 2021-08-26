#pragma once

#include <hydra/common.h>
#include <hydra/indexing/lintable.h>

namespace hydra {

template <class bit_t=std_bit_t> class Spinhalf {
public:
  Spinhalf() = default;
  Spinhalf(int n_sites, int n_up);

  inline int n_sites() const { return n_sites_; }

  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const {return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  
  inline idx_t size() const { return size_; }

  inline idx_t index(bit_t spins) const {
    return lintable_.index(spins);
  }

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

private:
  int n_sites_;

  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;

  LinTable<bit_t> lintable_;
  idx_t size_;
};

} // namespace hydra
