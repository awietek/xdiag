#include "spinhalf.h"

#include <hydra/combinatorics/binomial.h>
#include <hydra/combinatorics/combinations.h>
#include <hydra/combinatorics/subsets.h>

namespace hydra {
namespace detail {

int sz_to_nupspins(int sz, int n_sites) { return n_sites / 2 + sz; }

template <class bit_t, class BitsGen, class SymmetryGroup>
void fill_states_norms(BitsGen &&bits_generator, SymmetryGroup &&symmetry_group,
                       Representation const &irrep, std::vector<bit_t> &states,
                       std::vector<double> &norms) {
  for (auto bits : bits_generator) {
    bit_t rep = symmetry_group.representative(bits);

    if (rep == bits) {
      // Compute norm
      complex amplitude = 0.0;
      for (int sym = 0; sym < (int)irrep.characters.size(); ++sym) {
        bit_t derivedstate = symmetry_group.apply(sym, bits);
        if (derivedstate == bits)
          amplitude += irrep.characters[sym];
      }
      if (std::abs(amplitude) > 1e-8) {
        states.push_back(rep);
        norms.push_back(std::sqrt(std::abs(amplitude)));
      }
    }
  }
} // fill_states_norms
} // namespace detail

template <class bit_t, class SymmetryGroup>
Spinhalf<bit_t, SymmetryGroup>::Spinhalf(int n_sites)
    : n_sites_(n_sites), sz_conserved_(false), symmetry_group_defined_(false),
      size_((idx_t)1 << n_sites) {
  auto advance = [this](bit_t &state, idx_t &idx) {
    ++state;
    ++idx;
  };
  begin_ = SpinhalfIterator<bit_t>(0, 0, advance);
  end_ = SpinhalfIterator<bit_t>((bit_t)1 << n_sites, size_, advance);
}

template <class bit_t, class SymmetryGroup>
Spinhalf<bit_t, SymmetryGroup>::Spinhalf(int n_sites, int sz)
    : n_sites_(n_sites), sz_conserved_(true), sz_(sz),
      n_upspins_(detail::sz_to_nupspins(sz, n_sites)),
      lintable_(n_sites, n_upspins_), symmetry_group_defined_(false),
      size_(combinatorics::binomial(n_sites, n_upspins_)) {

  if ((n_upspins_ < 0) || (n_upspins_ > n_sites))
    HydraLog.err("Error creating Spinhalf: "
                 "invalid value of sz");

  auto advance = [this](bit_t &state, idx_t &idx) {
    state = combinatorics::get_next_pattern(state);
    ++idx;
  };
  bit_t b = ((bit_t)1 << n_upspins_) - 1;
  bit_t e = combinatorics::get_next_pattern(b << (n_sites_ - n_upspins_));
  begin_ = SpinhalfIterator<bit_t>(b, 0, advance);
  end_ = SpinhalfIterator<bit_t>(e, size_, advance);
}

template <class bit_t, class SymmetryGroup>
Spinhalf<bit_t, SymmetryGroup>::Spinhalf(int n_sites,
                                         SymmetryGroup symmetry_group,
                                         Representation irrep)
    : n_sites_(n_sites), sz_conserved_(false), symmetry_group_defined_(true),
      symmetry_group_(symmetry_group), irrep_(irrep) {

  detail::fill_states_norms(Subsets<bit_t>(n_sites), symmetry_group, irrep,
                            states_, norms_);
  size_ = (idx_t)states_.size();
  assert((idx_t)norms_.size() == size_);
  auto advance = [this](bit_t &state, idx_t &idx) {
    state = this->states_[++idx];
  };
  bit_t b = (size_ > 0) ? states_[0] : 0;
  bit_t e = (size_ > 0) ? states_[size_ - 1] : 0;
  begin_ = SpinhalfIterator<bit_t>(b, 0, advance);
  end_ = SpinhalfIterator<bit_t>(e, size_, advance);
}

template <class bit_t, class SymmetryGroup>
Spinhalf<bit_t, SymmetryGroup>::Spinhalf(int n_sites, int sz,
                                         SymmetryGroup symmetry_group,
                                         Representation irrep)
    : n_sites_(n_sites), sz_conserved_(true), sz_(sz),
      n_upspins_(detail::sz_to_nupspins(sz, n_sites)),
      symmetry_group_defined_(true), symmetry_group_(symmetry_group),
      irrep_(irrep) {
  if ((n_upspins_ < 0) || (n_upspins_ > n_sites))
    HydraLog.err("Error creating Spinhalf: "
                 "invalid value of sz");

  detail::fill_states_norms(Combinations<bit_t>(n_sites, n_upspins_),
                            symmetry_group, irrep, states_, norms_);

  size_ = (idx_t)states_.size();
  assert((idx_t)norms_.size() == size_);
  auto advance = [this](bit_t &state, idx_t &idx) {
    state = this->states_[++idx];
  };
  bit_t b = (size_ > 0) ? states_[0] : 0;
  bit_t e = (size_ > 0) ? states_[size_ - 1] : 0;
  begin_ = SpinhalfIterator<bit_t>(b, 0, advance);
  end_ = SpinhalfIterator<bit_t>(e, size_, advance);
}

template class Spinhalf<uint16, SpaceGroup<uint16>>;
template class Spinhalf<uint32, SpaceGroup<uint32>>;
template class Spinhalf<uint64, SpaceGroup<uint64>>;

} // namespace hydra
