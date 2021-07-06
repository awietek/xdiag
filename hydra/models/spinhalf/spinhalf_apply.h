#pragma once

#include <lila/all.h>

#include <hydra/common.h>
#include <hydra/models/spinhalf/spinhalf.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra {

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, lila::Vector<double> const &vec_in,
           Spinhalf<bit_t> const &block_out, lila::Vector<double> &vec_out);

template <class bit_t>
void Apply(BondList const &bonds, Couplings const &couplings,
           Spinhalf<bit_t> const &block_in, lila::Vector<complex> const &vec_in,
           Spinhalf<bit_t> const &block_out, lila::Vector<complex> &vec_out);

} // namespace hydra
