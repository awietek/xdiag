#include <iostream>
#include <hydra/utils/bitops.h>
#include <hydra/utils/complex.h>
#include "hubbardmodeldetail.h"
#include "hubbardmodelmpi.h">

hs_holes_ = Spinhalf<state_t>(n_sites_, n_sites_ - qn_.n_upspins - qn_.n_downspins);

indexing_holes_ = IndexTable<Spinhalf<state_t>, uint64>(hs_holes_);

inline uint64 hole_idx(state_t upspins, state_t downspins)
{ 
  return indexing_holes_.index(upspins & downspins);
}
