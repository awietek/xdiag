#pragma once

#include <hydra/common.h>
#include <hydra/parameters/parrameters.h>
#include <hydra/symmetries/charactertable.h>

namespace hydra {

  template <class bit_t>
class Spinhalf {
public:
  Spinhalf() = default;
  Spinhalf(int n_sites);
  Spinhalf(int n_sites, Parameters parameters);
  Spinhalf(int n_sites, CharacterTable character_table, Parameters parameters);

  int n_sites() const { return n_sites_; }
  idx_t dim() const { return dim_; }

private:
  int n_sites_;
  Parameters parameters_;
  CharacterTable character_table_;

    LinTable lintable_;

    std::vector<bit_t> states_;
    std::vector<double> norms_;
  };

lila::Matrix<double> matrix(BondList bondlist, Couplings couplings,
                            Spinhalf const &model);

void apply(BondList bondlist, Couplings couplings, Spinhalf const &block_in,
           lila::Vector<coeff_t> const &vec_in, Spinhalf const &block_out,
           lila::Vector<coeff_t> &vec_out)

} // namespace hydra
