#pragma once

#include <hydra/bitops/bitops.h>
#include <hydra/common.h>
#include <hydra/operators/bond.h>

namespace hydra::tj {

template <typename bit_t, typename coeff_t, class Indexing, class Filler>
void apply_symmetric_ising(Bond const &bond, Indexing &&indexing,
                           Filler &&fill) {
  assert(bond.coupling_defined());
  assert(bond.type_defined());
  assert((bond.type() == "ISING") || (bond.type() == "TJISING"));
  assert(bond.size() == 2);
  assert(bond.sites_disjoint());

  int n_sites = indexing.n_sites();

  coeff_t J = bond.coupling<coeff_t>();
  int s1 = bond[0];
  int s2 = bond[1];

  // Set values for same/diff (tJ model definition)
  coeff_t val_same, val_diff;
  if (bond.type() == "ISING") {
    val_same = J / 4.;
    val_diff = -J / 4.;
  } else { // (bond.type() == "TJISING")
    val_same = 0.;
    val_diff = -J / 2.;
  }

  // bitmasks for fast evaluations
  bit_t s1mask = (bit_t)1 << s1;
  bit_t s2mask = (bit_t)1 << s2;
  bit_t mask = s1mask | s2mask;
  bit_t sitesmask = ((bit_t)1 << n_sites) - 1;

  // Loop over all up configurations
  for (idx_t idx_up = 0; idx_up < indexing.n_rep_ups(); ++idx_up) {
    bit_t ups = indexing.rep_ups(idx_up);
    idx_t ups_offset = indexing.ups_offset(idx_up);
    auto dnss = indexing.dns_for_ups_rep(ups);

    if (bitops::popcnt(ups & mask) == 2) { // both spins pointing up
      for (idx_t idx = ups_offset; idx < ups_offset + (idx_t)dnss.size();
           ++idx) {
        fill(idx, idx, val_same);
      }
    } else {
      auto syms = indexing.syms_ups(ups);
      idx_t idx = ups_offset;

      // ups have trivial stabilizer => dns need to be deposited
      if (syms.size() == 1) {
        bit_t not_ups = (~ups) & sitesmask;

        if (ups & s1mask) { // s1 is pointing up

          for (auto dnsc : dnss) {
            bit_t dns = bitops::deposit(dnsc, not_ups);
            if (dns & s2mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else if (ups & s2mask) { // s2 is pointing up

          for (auto dnsc : dnss) {
            bit_t dns = bitops::deposit(dnsc, not_ups);
            if (dns & s1mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else { // no upspins

          for (auto dnsc : dnss) {
            bit_t dns = bitops::deposit(dnsc, not_ups);
            if (bitops::popcnt(dns & mask) == 2) {
              fill(idx, idx, val_same);
            }
            ++idx;
          }
        }
      }      // if (syms_up.size() == 1)
      else { // ups have stabilizer => dns don't need to be deposited
        if (ups & s1mask) { // s1 is pointing up

          for (auto dns : dnss) {
            if (dns & s2mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else if (ups & s2mask) { // s2 is pointing up

          for (auto dns : dnss) {
            if (dns & s1mask) {
              fill(idx, idx, val_diff);
            }
            ++idx;
          }

        } else { // no upspins

          for (auto dns : dnss) {
            if (bitops::popcnt(dns & mask) == 2) {
              fill(idx, idx, val_same);
            }
            ++idx;
          }
        }
      }

    } // else of if (popcnt(ups & mask) == 2)
  }   // loop over ups
} // loop over bonds

} // namespace hydra::tj
