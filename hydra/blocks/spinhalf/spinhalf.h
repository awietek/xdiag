#pragma once

#include <memory>

#include <hydra/common.h>

#include <hydra/indexing/indexing_variants.h>

#include <hydra/symmetries/permutation_group.h>
#include <hydra/symmetries/representation.h>

namespace hydra {

class Spinhalf {
public:
  using indexing_t = indexing::SpinhalfIndexing;

  Spinhalf() = default;
  Spinhalf(int64_t n_sites);
  Spinhalf(int64_t n_sites, int64_t n_up);
  Spinhalf(int64_t n_sites, PermutationGroup permutation_group,
           Representation irrep, int64_t n_sublat = 0);
  Spinhalf(int64_t n_sites, int64_t n_up, PermutationGroup permutation_group,
           Representation irrep, int64_t n_sublat = 0);

  inline int64_t n_sites() const { return n_sites_; }
  inline bool sz_conserved() const { return sz_conserved_; }
  inline int64_t sz() const { return sz_; }
  inline int64_t n_up() const { return n_up_; }
  inline int64_t n_dn() const { return n_dn_; }
  inline bool symmetric() const { return symmetric_; }
  inline PermutationGroup const &permutation_group() const {
    return permutation_group_;
  }
  inline Representation const &irrep() const { return irrep_; }
  inline idx_t size() const { return size_; }

  bool operator==(Spinhalf const &rhs) const;
  bool operator!=(Spinhalf const &rhs) const;

  indexing_t const &indexing() const;

private:
  int64_t n_sites_;
  bool sz_conserved_;
  int64_t n_up_;
  int64_t n_dn_;
  int64_t sz_;
  bool symmetric_;
  int64_t n_sublat_;
  PermutationGroup permutation_group_;
  Representation irrep_;
  std::shared_ptr<indexing_t> indexing_;
  idx_t size_;
};

} // namespace hydra
