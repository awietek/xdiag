#pragma once
#include "extern/armadillo/armadillo"

#include <hydra/common.h>
#include <hydra/indexing/indexing_variants.h>
#include <hydra/operators/bondlist.h>
#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

class tJ {
public:
  using indexing_t = indexing::tJIndexing;

  tJ() = default;
  tJ(int64_t n_sites, int64_t nup, int64_t ndn);
  tJ(int64_t n_sites, int64_t nup, int64_t ndn, PermutationGroup permutation_group,
     Representation irrep);

  inline int64_t n_sites() const { return n_sites_; }
  inline int64_t n_up() const { return n_up_; }
  inline int64_t n_dn() const { return n_dn_; }

  inline bool charge_conserved() const { return charge_conserved_; }
  inline bool sz_conserved() const { return sz_conserved_; }

  inline bool symmetric() const { return symmetric_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  bool operator==(tJ const &rhs) const;
  bool operator!=(tJ const &rhs) const;

  indexing_t const &indexing() const;

private:
  int64_t n_sites_;
  bool charge_conserved_;
  int64_t charge_;
  bool sz_conserved_;
  int64_t sz_;
  int64_t n_up_;
  int64_t n_dn_;
  bool symmetric_;
  PermutationGroup permutation_group_;
  Representation irrep_;
  std::shared_ptr<indexing_t> indexing_;
  idx_t size_;
};

} // namespace hydra
