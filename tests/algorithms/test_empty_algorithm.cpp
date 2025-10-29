// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "../catch.hpp"

#include <iostream>

#include "../blocks/electron/testcases_electron.hpp"
#include "../blocks/spinhalf/testcases_spinhalf.hpp"

#include <xdiag/algorithms/lanczos/eigs_lanczos.hpp>
#include <xdiag/algorithms/lanczos/eigvals_lanczos.hpp>
#include <xdiag/algorithms/sparse_diag.hpp>
#include <xdiag/algorithms/time_evolution/imaginary_time_evolve.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve.hpp>
#include <xdiag/algorithms/time_evolution/evolve_lanczos.hpp>
#include <xdiag/algorithms/time_evolution/time_evolve_expokit.hpp>

#include <xdiag/states/create_state.hpp>
#include <xdiag/utils/xdiag_show.hpp>

using namespace xdiag;

TEST_CASE("empty_algorithm", "[algorithm]") try {
  using namespace testcases::spinhalf;
  using testcases::electron::get_cyclic_group_irreps;
  Log("empty algorithms");
  auto irreps = get_cyclic_group_irreps(4);
  auto ops = HBchain(4, 1.0, 1.0);

  // create empty block
  auto block = Spinhalf(4, 4, irreps[2]);

  auto rstate = random_state(block);
  double e0 = eigval0(ops, block);
  auto [e00, psi0] = eig0(ops, block);

  auto res1 = eigvals_lanczos(ops, block);
  auto res2 = eigvals_lanczos(ops, rstate);
  auto res3 = eigvals_lanczos_inplace(ops, rstate);

  auto rres1 = eigs_lanczos(ops, block);
  auto rres2 = eigs_lanczos(ops, rstate);

  auto psi1 = time_evolve(ops, rstate, 1.0);
  time_evolve_inplace(ops, rstate, 1.0);

  auto psi2 = imaginary_time_evolve(ops, rstate, 1.0);
  imaginary_time_evolve_inplace(ops, rstate, 1.0);

  auto rrres1 = evolve_lanczos(ops, rstate, 1.0);
  auto rrres2 = evolve_lanczos(ops, rstate, complex(1.0, 1.0));
  auto rrres3 = evolve_lanczos_inplace(ops, rstate, 1.0);
  auto rrres4 = evolve_lanczos_inplace(ops, rstate, complex(1.0, 1.0));


  auto rrres5 = time_evolve_expokit(ops, rstate, 1.0);
  auto rrres6 = time_evolve_expokit_inplace(ops, rstate, 1.0);

} catch (Error e) {
  error_trace(e);
}
