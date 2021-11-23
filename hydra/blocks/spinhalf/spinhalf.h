#pragma once

#include <hydra/common.h>

#include <hydra/blocks/spinhalf/spinhalf_apply.h>
#include <hydra/blocks/spinhalf/spinhalf_matrix.h>

#include <hydra/indexing/spinhalf/spinhalf_indexing.h>

#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <typename bit_t> class Spinhalf {
public:
  Spinhalf() = default;
  Spinhalf(int n_sites, int n_up);

  inline int n_sites() const { return n_sites_; }
  inline bool sz_conserved() const { return sz_conserved_; }
  inline int sz() const { return sz_; }
  inline int n_up() const { return n_up_; }
  inline int n_dn() const { return n_dn_; }
  inline idx_t size() const { return size_; }

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

  indexing::SpinhalfIndexing<bit_t> const &indexing() const {
    return indexing_;
  }

private:
  int n_sites_;
  bool sz_conserved_;
  int n_up_;
  int n_dn_;
  int sz_;
  indexing::SpinhalfIndexing<bit_t> indexing_;
  idx_t size_;
};

} // namespace hydra
