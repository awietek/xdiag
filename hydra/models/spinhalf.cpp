#include "spinhalf.h"

namespace hydra {

template <class bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites)
    : n_sites_(n_sites), sz_conserved_(false), sz_(0),
      spinflip_symmetric_(false), spinflip_sign_(0.0),
      space_group_defined_(false) {}

template <class bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, Parameters p)
    : n_sites_(n_sites), sz_conserved_(p.defined("Sz")),
      sz_(p.defined("Sz") ? (int)p["Sz"] : 0),
      spinflip_symmetric_(p.defined("si")),
      spinflip_sign_(p.defined("si") ? (int)p["si"] : 0),
      space_group_defined_(false), lintable_(n_sites) {}

namespace detail {

template <class bit_t, class BitsGen, class FindRep, class ApplySym>
void fill_states_norms(BitsGen &&bits_generator, FindRep &&find_rep,
                       ApplySym &&apply_sym, std::vector<complex> characters,
                       std::vector<bit_t> &states, std::vector<bit_t> &norms) {
  for (auto bits : bits_generator) {
    bit_t rep = find_rep(bits);
    if (rep == bits) {
      // Compute norm
      complex amplitude = 0.0;
      for (int sym = 0; sym < (int)characters.size(); ++sym) {
        bit_t derivedstate = apply_sym(bits, sym);
        if (derivedstate == bits)
          amplitude += characters[sym];
      }

      if (std::abs(amplitude) < 1e-8) {
        states.push_back(rep);
        norms.push_back(std::sqrt(std::abs(amplitude)));
      }

    } // (rep == bits)
  }   // for (auto bits : bits_generator)
} // fill_states_norms
} // namespace detail

template <class bit_t>
Spinhalf<bit_t>::Spinhalf(int n_sites, CharacterTable character_table,
                          Parameters p)
    : n_sites_(n_sites), sz_conserved_(p.defined("Sz")),
      sz_(p.defined("Sz") ? (int)p["Sz"] : 0),
      spinflip_symmetric_(p.defined("si")),
      spinflip_sign_(p.defined("si") ? (int)p["si"] : 0),
      space_group_defined_(p.defined("k")),
      space_group_irrep_(p.defined("k") ? (std::string)p["k"] : ""),
      character_table_(character_table),
      space_group_(character_table_.space_group()) {
  assert(n_sites >= 0);

  // Create list of representatives if spacegroup is defined
  if (space_group_defined_) {

    if (spinflip_symmetric_) {
      assert((spinflip_sign_ == 1) || (spinflip_sign_ == -1));

      // Prepare the characters for the symmetry
      int n_space_sym = space_group_.n_symmetries();
      double si = (double)spinflip_sign_;
      for (int sym = 0; sym < n_space_sym; ++sym) {
        complex sg_char = character_table_(space_group_irrep_, sym);
        characters_.push_back(sg_char);
        characters_.push_back(sg_char * si);
      }

      bit_t mask = ((bit_t)1 << n_sites) - 1;

      // Define functions find representative and apply symmetries
      auto find_rep = [&mask, this](bit_t bits) -> bit_t {
        bit_t rep = this->space_group_.representative(bits);
        return std::min(rep, (~rep) & mask);
      };
      auto apply_sym = [&mask, this](bit_t bits, int sym) -> bit_t {
        bit_t dbits = this->space_group_.apply_symmetry(bits, sym >> 1);
        if (bits & 1)
          return dbits;
        else
          return (~dbits) & mask;
      };

      if (sz_conserved_)
        detail::fill_states_norms(Combinations(n_sites_, sz_), find_rep,
                                  apply_sym, characters_, states_, norms_);
      else
        detail::fill_states_norms(Subsets(n_sites_), find_rep, apply_sym,
                                  characters_, states_, norms_);

    } else { // not spinflip_symmetric_
      characters_ = character_table.characters(space_group_irrep_);

      // Define functions find representative and apply symmetries
      auto find_rep = [this](bit_t bits) -> bit_t {
        return this->space_group_.representative(bits);
      };
      auto apply_sym = [this](bit_t bits, int sym) -> bit_t {
        return this->space_group_.apply_symmetry(bits, sym);
      };

      if (sz_conserved_)
        detail::fill_states_norms(Combinations(n_sites_, sz_), find_rep,
                                  apply_sym, characters_, states_, norms_);
      else
        detail::fill_states_norms(Subsets(n_sites_), find_rep, apply_sym,
                                  characters_, states_, norms_);
    }

    size_ = (idx_t)states_.size();
    auto advance = [this](bit_t& state, idx_t& idx){
		     state = this->states[++idx];
		   };
    
    begin_ = SpinhalfIterator(states_[0], 0, advance);
    end_ = SpinhalfIterator(states_[size_-1], size_, advance);
  } else { // not space_group_defined_
  }
}

} // namespace hydra
