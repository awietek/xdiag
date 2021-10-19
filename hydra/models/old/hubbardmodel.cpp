#include <lila/all.h>

#include <hydra/bases/basis_spinhalf.h>
#include <hydra/indexing/index_electron.h>
#include <hydra/indexing/index_table.h>
#include <hydra/utils/bitops.h>
#include <hydra/utils/range.h>

#include "hubbardmodel.h"
#include "hubbardmodeldetail.h"

namespace hydra {

template <class coeff_t, class bit_t, class idx_t>
HubbardModel<coeff_t, bit_t, idx_t>::HubbardModel(BondList bondlist,
                                                  Couplings couplings, qn_t qn)
    : n_sites_(bondlist.n_sites()), qn_(qn) {
  BasisElectron<bit_t> hs(n_sites_, qn_);
  dim_ = hs.size();

  detail::set_hubbard_terms<coeff_t>(
      bondlist, couplings, hoppings_, hopping_amplitudes_, currents_,
      current_amplitudes_, interactions_, interaction_strengths_, onsites_,
      onsite_potentials_, szszs_, szsz_amplitudes_, exchanges_,
      exchange_amplitudes_, U_);
}

template <class coeff_t, class bit_t, class idx_t>
void HubbardModel<coeff_t, bit_t, idx_t>::set_qn(qn_t qn) {
  qn_ = qn;
  BasisElectron<bit_t> hs(n_sites_, qn_);
  dim_ = hs.size();
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<coeff_t> HubbardModel<coeff_t, bit_t, idx_t>::matrix() const {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  BasisElectron<bit_t> hs(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing(hs);
  int dim = indexing.size();
  lila::Matrix<coeff_t> hamilton(dim, dim);
  lila::Zeros(hamilton);

  // Apply Hubbard U term
  if (std::abs(U_) > 1e-14) {
    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      hamilton(idx, idx) += U_ * (double)popcnt(state.ups & state.dns);
    }
  }

  // Apply Hubbard V term
  int interaction_idx = 0;
  for (auto pair : interactions_) {
    const int s1 = pair.first;
    const int s2 = pair.second;
    const double V = interaction_strengths_[interaction_idx];
    if (std::abs(V) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        hamilton(idx, idx) +=
            V * (double)((gbit(state.ups, s1) + gbit(state.dns, s1)) *
                         (gbit(state.ups, s2) + gbit(state.dns, s2)));
      }
    }
    ++interaction_idx;
  }

  // Apply onsite chemical potential
  int onsite_idx = 0;
  for (auto site : onsites_) {
    const double mu = onsite_potentials_[onsite_idx];
    if (std::abs(mu) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        hamilton(idx, idx) -=
            mu * (double)((gbit(state.ups, site) + gbit(state.dns, site)));
      }
    }
    ++onsite_idx;
  }

  // Apply hopping terms
  int hopping_idx = 0;
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t t = hopping_amplitudes_[hopping_idx];
    if (pair.first > pair.second)
      t = lila::conj(t);

    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    bit_t secondmask = (bit_t)1 << s2;

    if (std::abs(t) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        const bit_t &ups = state.ups;
        const bit_t &dns = state.dns;

        // upspins hopping
        if (((ups & flipmask) != 0) && ((ups & flipmask) != flipmask)) {
          const double fermi =
              popcnt(gbits(ups ^ dns, s2 - s1, s1)) & 1 ? 1. : -1.;
          auto new_state = state;
          new_state.ups = ups ^ flipmask;
          int new_idx = indexing.index(new_state);
          if (ups & secondmask)
            hamilton(new_idx, idx) += fermi * t;
          else
            hamilton(new_idx, idx) -= lila::conj(t) * fermi;
        }

        // dns hopping
        if (((dns & flipmask) != 0) && ((dns & flipmask) != flipmask)) {
          const double fermi =
              popcnt(gbits(ups ^ (dns ^ secondmask), s2 - s1, s1 + 1)) & 1
                  ? 1.
                  : -1.;
          auto new_state = state;
          new_state.dns = dns ^ flipmask;
          int new_idx = indexing.index(new_state);
          if (dns & secondmask)
            hamilton(new_idx, idx) += fermi * t;
          else
            hamilton(new_idx, idx) -= lila::conj(t) * fermi;
        }
      }
    }
    ++hopping_idx;
  } // loop over hoppings

  // SzSz terms
  int szsz_idx = 0;
  for (auto pair : szszs_) {
    int s1 = pair.first;
    int s2 = pair.second;
    double jz = szsz_amplitudes_[szsz_idx] * 0.25;

    if (std::abs(jz) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        const bit_t &ups = state.ups;
        const bit_t &dns = state.dns;

        bool up1 = gbit(ups, s1);
        bool up2 = gbit(ups, s2);
        bool dn1 = gbit(dns, s1);
        bool dn2 = gbit(dns, s2);
        //      			auto coeff = jz*(((double)gbit(ups,
        //      s1) - (double)gbit(dns, s1)) *
        //			    ((double)gbit(ups, s2) -
        //(double)gbit(dns, s2)));

        if (!(up1 && dn1) && !(dn2 && up2)) {
          if ((up1 && up2) || (dn1 && dn2))
            hamilton(idx, idx) += jz;
          else if ((up1 && dn2) || (dn1 && up2))
            hamilton(idx, idx) += -jz;
        }
      }
    }
    ++szsz_idx;
  }

  // Exchange terms
  int exchange_idx = 0;
  for (auto pair : exchanges_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t jx = exchange_amplitudes_[exchange_idx] * 0.5;
    const bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    if (std::abs(jx) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        const bit_t &ups = state.ups;
        const bit_t &dns = state.dns;

        if ((popcnt(ups & flipmask) == 1) && (popcnt(dns & flipmask) == 1) &&
            popcnt((dns & flipmask) & (ups & flipmask)) == 0)
        // if ((up1 && dn2 && !up2 && !dn1) || (dn1 && up2 && !dn2 && !up1))
        {
          auto new_state = state;
          new_state.ups = ups ^ flipmask;
          new_state.dns = dns ^ flipmask;

          // if (popcnt(new_state.ups) != popcnt(state.ups)) {
          //   std::cout << PrintSpinhalf(32, state.ups) << std::endl
          //             << s1 << " " << s2 << std::endl
          //             << popcnt(state.ups) << std::endl
          //             << PrintSpinhalf(32, new_state.ups) << std::endl
          //             << popcnt(new_state.ups) << std::endl
          //             << std::endl;
          // }
          // assert(popcnt(new_state.ups) == popcnt(state.ups));
          // assert(popcnt(new_state.dns) == popcnt(state.dns));

          int flipped_idx = indexing.index(new_state);
          hamilton(flipped_idx, idx) += jx;
        }

      } // loop over spin configurations
    }   // if (std::abs(jx) > 1e-14)
    ++exchange_idx;
  } // for (auto pair: exchanges_)

  return hamilton;
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<double>
HubbardModel<coeff_t, bit_t, idx_t>::szMatrix(int siteIndex) const {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  assert(0 <= siteIndex && siteIndex < n_sites_);

  BasisElectron<bit_t> hs(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing(hs);
  int dim = indexing.size();
  lila::Matrix<double> sz(dim, dim);
  lila::Zeros(sz);

  // Assemble sz matrix

  for (int idx : range<>(indexing.size())) {
    auto state = indexing.state(idx);
    const bit_t &ups = state.ups;
    const bit_t &dns = state.dns;
    bool up = gbit(ups, siteIndex);
    bool dn = gbit(dns, siteIndex);

    sz(idx, idx) = up - dn;
  }
  return sz;
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<double>
HubbardModel<coeff_t, bit_t, idx_t>::sPlusMatrix(int siteIndex) const {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  assert(0 <= siteIndex && siteIndex < n_sites_);

  BasisElectron<bit_t> hs(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing(hs);

  if (qn_.n_up < n_sites_ && qn_.n_dn > 0) {

    // Initialize Hilbert space for target space
    auto target_qn = qn_;
    ++target_qn.n_up;
    --target_qn.n_dn;
    BasisElectron<bit_t> target_hs(n_sites_, target_qn);
    IndexElectron<bit_t, idx_t> target_indexing(target_hs);

    int base_dim = indexing.size();
    int target_dim = target_indexing.size();
    lila::Matrix<double> sPlus;
    sPlus.resize(target_dim, base_dim);
    lila::Zeros(sPlus);

    // Assemble sPlus matrix

    bit_t flipmask = ((bit_t)1 << siteIndex);

    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      const bit_t &ups = state.ups;
      const bit_t &dns = state.dns;
      bool up = gbit(ups, siteIndex);
      bool dn = gbit(dns, siteIndex);

      if (dn && !up) {
        auto flipped_state = state;
        flipped_state.ups = state.ups ^ flipmask;
        flipped_state.dns = state.dns ^ flipmask;
        int target_idx = target_indexing.index(flipped_state);

        sPlus(target_idx, idx) = 1;
      }
    }
    return sPlus;
  } else {
    lila::Matrix<double> sPlus(1, 1);
    lila::Zeros(sPlus);
    return sPlus;
  }
}

template <class coeff_t, class bit_t, class idx_t>
lila::Matrix<double>
HubbardModel<coeff_t, bit_t, idx_t>::sMinusMatrix(int siteIndex) const {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  assert(0 <= siteIndex && siteIndex < n_sites_);

  BasisElectron<bit_t> hs(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing(hs);

  if (qn_.n_dn < n_sites_ && qn_.n_up > 0) {

    // Initialize Hilbert space for target space
    auto target_qn = qn_;
    --target_qn.n_up;
    ++target_qn.n_dn;
    BasisElectron<bit_t> target_hs(n_sites_, target_qn);
    IndexElectron<bit_t, idx_t> target_indexing(target_hs);

    int base_dim = indexing.size();
    int target_dim = target_indexing.size();
    lila::Matrix<double> sMinus;
    sMinus.resize(target_dim, base_dim);
    lila::Zeros(sMinus);

    // Assemble sMinus matrix

    bit_t flipmask = ((bit_t)1 << siteIndex);

    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      const bit_t &ups = state.ups;
      const bit_t &dns = state.dns;
      bool up = gbit(ups, siteIndex);
      bool dn = gbit(dns, siteIndex);

      if (!dn && up) {

        auto flipped_state = state;
        flipped_state.ups = state.ups ^ flipmask;
        flipped_state.dns = state.dns ^ flipmask;
        int target_idx = target_indexing.index(flipped_state);

        sMinus(target_idx, idx) = 1;
      }
    }
    return sMinus;
  } else {
    lila::Matrix<double> sMinus(1, 1);
    lila::Zeros(sMinus);
    return sMinus;
  }
}

template <class coeff_t, class bit_t, class idx_t>
void HubbardModel<coeff_t, bit_t, idx_t>::apply_hamiltonian(
    const lila::Vector<coeff_t> &in_vec, lila::Vector<coeff_t> &out_vec) const {

  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  BasisElectron<bit_t> hs(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing(hs);
  idx_t dim = indexing.size();
  assert(in_vec.size() == dim);
  assert(out_vec.size() == dim);
  lila::Zeros(out_vec);

  // Apply Hubbard U term
  if (std::abs(U_) > 1e-14) {
    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      auto coeff = U_ * (double)popcnt(state.ups & state.dns);
      out_vec(idx) += coeff * in_vec(idx);
    }
  }

  // Apply Hubbard V term
  int interaction_idx = 0;
  for (auto pair : interactions_) {
    const int s1 = pair.first;
    const int s2 = pair.second;
    const double V = interaction_strengths_[interaction_idx];
    if (std::abs(V) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        auto coeff = V * (double)((gbit(state.ups, s1) + gbit(state.dns, s1)) *
                                  (gbit(state.ups, s2) + gbit(state.dns, s2)));
        out_vec(idx) += coeff * in_vec(idx);
      }
    }
    ++interaction_idx;
  }

  // Apply onsite chemical potential
  int onsite_idx = 0;
  for (auto site : onsites_) {
    const double mu = onsite_potentials_[onsite_idx];
    if (std::abs(mu) > 1e-14) {
      for (int idx : range<>(indexing.size())) {
        auto state = indexing.state(idx);
        auto coeff =
            mu * (double)((gbit(state.ups, site) + gbit(state.dns, site)));
        out_vec(idx) -= coeff * in_vec(idx);
      }
    }
    ++onsite_idx;
  }

  // Apply hopping terms
  int hopping_idx = 0;
  for (auto pair : hoppings_) {
    int s1 = std::min(pair.first, pair.second);
    int s2 = std::max(pair.first, pair.second);
    coeff_t t = hopping_amplitudes_[hopping_idx];
    if (pair.first > pair.second)
      t = lila::conj(t);

    bit_t flipmask = ((bit_t)1 << s1) | ((bit_t)1 << s2);
    for (int idx : range<>(indexing.size())) {
      auto state = indexing.state(idx);
      const bit_t &ups = state.ups;
      const bit_t &dns = state.dns;

      // ups hopping
      if (((ups & flipmask) != 0) && ((ups & flipmask) != flipmask)) {
        const double fermi =
            popcnt(gbits(ups, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1. : -1.;
        auto new_state = state;
        new_state.ups = ups ^ flipmask;
        int new_idx = indexing.index(new_state);
        out_vec(new_idx) -= fermi * t * in_vec(idx);
      }

      // dns hopping
      if (((dns & flipmask) != 0) && ((dns & flipmask) != flipmask)) {
        const double fermi =
            popcnt(gbits(dns, s2 - s1 - 1, s1 + 1)) % 2 == 0 ? 1. : -1.;
        auto new_state = state;
        new_state.dns = dns ^ flipmask;
        int new_idx = indexing.index(new_state);
        out_vec(new_idx) -= fermi * t * in_vec(idx);
      }
    }
    ++hopping_idx;
  } // loop over hoppings
}

template <class coeff_t, class bit_t, class idx_t>
qn_electron HubbardModel<coeff_t, bit_t, idx_t>::apply_fermion(
    const lila::Vector<coeff_t> &state_before,
    lila::Vector<coeff_t> &state_after, std::string type, int site) const {
  using utils::gbit;
  using utils::gbits;
  using utils::popcnt;
  using utils::range;

  BasisElectron<bit_t> hs_before(n_sites_, qn_);
  IndexElectron<bit_t, idx_t> indexing_before(hs_before);

  auto qn_after = qn_;
  if (type == "cdagup")
    ++qn_after.n_up;
  else if (type == "cup")
    --qn_after.n_up;
  else if (type == "cdagdn")
    ++qn_after.n_dn;
  else if (type == "cdn")
    --qn_after.n_dn;
  else {
    std::cerr << "Error in apply_fermion: Invalid fermion type!" << std::endl;
    exit(EXIT_FAILURE);
  }

  BasisElectron<bit_t> hs_after(n_sites_, qn_after);
  IndexElectron<bit_t, idx_t> indexing_after(hs_after);
  int64 dim_after = indexing_after.size();

  state_after.resize(dim_after);

  const bit_t sitemask = (1 << site);
  const bit_t antisitemask = ~(1 << site);
  if (type == "cdagup") {
    for (int idx : range<>(indexing_before.size())) {
      auto config = indexing_before.state(idx);
      const bit_t &ups = config.ups;

      // raise local site val if 0
      if (gbit(ups, site) == 0) {
        auto new_config = config;
        new_config.ups |= sitemask;
        int new_idx = indexing_after.index(new_config);
        double fermi = popcnt(gbits(ups, site, 0)) % 2 == 0 ? 1. : -1.;
        state_after(new_idx) += fermi * state_before(idx);
      }
    }
  } else if (type == "cup") {
    for (int idx : range<>(indexing_before.size())) {
      auto config = indexing_before.state(idx);
      const bit_t &ups = config.ups;

      // lower local site val if 1
      if (gbit(ups, site) == 1) {
        auto new_config = config;
        new_config.ups &= antisitemask;
        int new_idx = indexing_after.index(new_config);
        double fermi = popcnt(gbits(ups, site, 0)) % 2 == 0 ? 1. : -1.;
        state_after(new_idx) += fermi * state_before(idx);
      }
    }
  } else if (type == "cdagdn") {
    for (int idx : range<>(indexing_before.size())) {
      auto config = indexing_before.state(idx);
      const bit_t &dns = config.dns;

      // raise local site val if 0
      if (gbit(dns, site) == 0) {
        auto new_config = config;
        new_config.dns |= sitemask;
        int new_idx = indexing_after.index(new_config);
        double fermi = popcnt(gbits(dns, site, 0)) % 2 == 0 ? 1. : -1.;
        state_after(new_idx) += fermi * state_before(idx);
      }
    }
  } else if (type == "cdn") {
    for (int idx : range<>(indexing_before.size())) {
      auto config = indexing_before.state(idx);
      const bit_t &dns = config.dns;

      // lower local site val if 1
      if (gbit(dns, site) == 1) {
        auto new_config = config;
        new_config.dns &= antisitemask;
        int new_idx = indexing_after.index(new_config);
        double fermi = popcnt(gbits(dns, site, 0)) % 2 == 0 ? 1. : -1.;
        state_after(new_idx) += fermi * state_before(idx);
      }
    }
  }
  return qn_after;
}

template class HubbardModel<double, uint32, uint16>;
template class HubbardModel<double, uint32, uint32>;
template class HubbardModel<double, uint32, uint64>;

template class HubbardModel<double, uint64, uint16>;
template class HubbardModel<double, uint64, uint32>;
template class HubbardModel<double, uint64, uint64>;

template class HubbardModel<complex, uint32, uint16>;
template class HubbardModel<complex, uint32, uint32>;
template class HubbardModel<complex, uint32, uint64>;

template class HubbardModel<complex, uint64, uint16>;
template class HubbardModel<complex, uint64, uint32>;
template class HubbardModel<complex, uint64, uint64>;

} // namespace hydra
