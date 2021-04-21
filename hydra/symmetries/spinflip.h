#pragma once

#include <vector>
#include <hydra/common.h>

namespace hydra {

template <class bit_t = std_bit_t>
class Spinflip {
public:
  Spinflip(int n_sites, bool identity_only=false);

  inline bit_t apply(int sym, bit_t state) const {
    return sym ? (~state) & sitemask_ : state;
  }
  inline bit_t representative(bit_t state) const {
    return std::min((~state) & sitemask_, state);
  }
  inline std::tuple<bit_t, int> representative_index(bit_t state) const {
    bit_t fstate = (~state) & sitemask_;
    if (fstate < state) return {fstate, 1};
    else return {state, 0};
  }

  Spinflip<bit_t> subgroup(std::vector<int> const &symmetry_numbers) const;

  inline int n_sites() const { return n_sites_; }
  inline int n_symmetries() const { return n_symmetries_; }
  inline int permutation(int sym, int site) const { return site; }

private:
  int n_sites_;
  bit_t sitemask_;
  int n_symmetries_;
  std::vector<int> permutation_array_; // size = n_symmetries_*n_sites_
};

} // namespace hydra
