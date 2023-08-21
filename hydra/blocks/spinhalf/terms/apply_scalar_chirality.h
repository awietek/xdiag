#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>

#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_no_sym.h>
#include <hydra/blocks/spinhalf/terms/apply_term_offdiag_sym.h>

namespace hydra::spinhalf {

// Scalar chirality term: J S*(S x S)

template <typename bit_t, typename coeff_t, bool symmetric, class IndexingIn,
          class IndexingOut, class Fill>
void apply_scalar_chirality(Bond const &bond, IndexingIn &&indexing_in,
                            IndexingOut &&indexing_out, Fill &&fill) {
  using bitops::gbit;

  assert(bond.coupling_defined());
  assert(bond.type_defined() && (bond.type() == "SCALARCHIRALITY"));
  assert(bond.size() == 3);
  assert(bond.sites_disjoint());
  assert(is_complex<coeff_t>());

  complex J = bond.coupling();
  coeff_t Jquarter = 0.;
  coeff_t Jquarter_conj = 0.;
  if constexpr (is_real<coeff_t>()) {
    Log.err("Error in spinhalf::apply_scalar_chirality: scalar chirality term "
            "cannot be used with real coefficients");
    (void)J;
  } else {
    Jquarter = complex(0, -0.25) * J;
    Jquarter_conj = hydra::conj(Jquarter);
  }

  int64_t s1 = bond[0];
  int64_t s2 = bond[1];
  int64_t s3 = bond[2];
  bit_t spinmask = ((bit_t)1 << s1) | ((bit_t)1 << s2) | ((bit_t)1 << s3);

  // scalar chirality annihilates 000 and 111
  auto non_zero_term = [&spinmask](bit_t spins) -> bool {
    bit_t threespins = spins & spinmask;
    return (threespins != 0) && (threespins != spinmask);
  };

  // rotate three sites cyclic
  auto term_action_cyclic = [&spinmask, &s1, &s2, &s3, &Jquarter](
                                bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t b1 = gbit(spins, s1);
    bit_t b2 = gbit(spins, s2);
    bit_t b3 = gbit(spins, s3);
    bit_t threespins_cyclic = (b1 << s2) | (b2 << s3) | (b3 << s1);
    bit_t spins_void = spins & (~spinmask);
    bit_t spins_cyclic = spins_void | threespins_cyclic;
    return {spins_cyclic, Jquarter};
  };

  // rotate three sites acyclic
  auto term_action_acyclic = [&spinmask, &s1, &s2, &s3, &Jquarter_conj](
                                 bit_t spins) -> std::pair<bit_t, coeff_t> {
    bit_t b1 = gbit(spins, s1);
    bit_t b2 = gbit(spins, s2);
    bit_t b3 = gbit(spins, s3);
    bit_t threespins_acyclic = (b1 << s3) | (b2 << s1) | (b3 << s2);
    bit_t spins_void = spins & (~spinmask);
    bit_t spins_acyclic = spins_void | threespins_acyclic;
    return {spins_acyclic, Jquarter_conj};
  };

  // Dispatch either symmetric of unsymmetric term application
  if constexpr (symmetric) {
    spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action_cyclic, fill);
    spinhalf::apply_term_offdiag_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action_acyclic, fill);

  } else {
    spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action_cyclic, fill);
    spinhalf::apply_term_offdiag_no_sym<bit_t, coeff_t>(
        indexing_in, indexing_out, non_zero_term, term_action_acyclic, fill);
  }
}

} // namespace hydra::spinhalf
