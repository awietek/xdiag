#ifndef HYDRA_INDEXING_INDEXSYMMETRIZED_
#define HYDRA_INDEXING_INDEXSYMMETRIZED_

#include <vector>

#include <hydra/common.h>
#include <hydra/symmetries/charactertable.h>

namespace hydra {

template <class basis_t, class idx_t = std_idx_t, bool fermionic = false>
class IndexSymmetrized {
public:
  using state_t = typename basis_t::state_t;

  IndexSymmetrized() = default;
  IndexSymmetrized(const basis_t &basis,
                   const CharacterTable &character_table,
                   std::string representation_name);

  idx_t index(const state_t &state) const;
  std::string representation_name() const { return representation_name_; }
  state_t state(const idx_t &index) const { return states_[index]; }
  double norm(const idx_t &index) const { return norms_[index]; }
  idx_t size() const { return (idx_t)states_.size(); }

  int find_representative(state_t *state) const;

private:
  CharacterTable character_table_;
  std::string representation_name_;
  int n_symmetries_;
  SpaceGroup space_group_;
  std::vector<complex> characters_;

  std::vector<state_t> states_;
  std::vector<double> norms_;
};

} // namespace hydra

#endif
