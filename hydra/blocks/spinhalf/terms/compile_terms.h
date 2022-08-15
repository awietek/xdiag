#pragma once

#include <utility>

#include <hydra/common.h>
#include <hydra/operators/bondlist.h>
#include <hydra/operators/couplings.h>

namespace hydra::terms::spinhalf {

std::pair<BondList, Couplings> compile_terms(BondList const &bonds,
                                             Couplings const &couplings);
} // namespace hydra::terms::spinhalf
