#include <cstdlib>
#include <itensor/all.h>
#include <xdiag/all.h>
#include <iomanip>

#include "../testcases_entanglement_entropy.h"

int main(int argc, char *argv[]) {
  using namespace itensor;
  using namespace xdiag;
  using namespace xdiag::entanglementtestcases;
  
  std::cout << std::setprecision(16);

  int N = 3;
  BondList bondlist;
  Couplings couplings;
  std::tie(bondlist, couplings) = tJrandom_model(N);
  
  auto sites = tJ(N);
  auto ampo = AutoMPO(sites);
  seedRNG(2);
  
  for (auto bond : bondlist)
    {
      auto type = bond.type();
      auto cpl_name = bond.coupling();
      auto cpl = couplings[cpl_name];

      //////////////////////////////
      cpl = std::real(cpl);
      //////////////////////////////
      
      auto s1 = bond.sites()[0] + 1;
      auto s2 = bond.sites()[1] + 1;
  
      if (bond.type() == "HUBBARDHOP")
	{
	  ampo += -cpl,"Cdagup", s1, "Cup", s2;
	  ampo += -std::conj(cpl),"Cdagup", s2, "Cup", s1;
	  ampo += -cpl,"Cdagdn", s1, "Cdn", s2;
	  ampo += -std::conj(cpl),"Cdagdn", s2, "Cdn", s1;
	  // std::cout << "Hop " << cpl << " " << std::conj(cpl) << std::endl;
	}

      if (bond.type() == "HEISENBERG")
	{
	  // std::cout << "Ex " << cpl << std::endl;
	  ampo += cpl    , "Sz", s1, "Sz", s2;
	  ampo += 0.5*cpl, "S+", s1, "S-", s2;
	  ampo += 0.5*cpl, "S-", s1, "S+", s2;
	}
    }
  auto H = toMPO(ampo);

  auto sweeps = Sweeps(20);
  sweeps.maxdim() = 10, 40, 100, 200, 200;
  sweeps.noise() = 10, 1, 1E-1, 1E-2, 1E-16;
  sweeps.cutoff() = 1E-12;

  seedRNG(1);
  
  for (int n_up = 0; n_up <= N; ++n_up)
    for (int n_dn = 0; n_dn <= N; ++n_dn) {

      // if (n_up + n_dn <= N) {
      if ((n_up ==1 && n_dn == 1)) {

        // Create a random starting state
        auto state = InitState(sites);
	for (auto i : range1(N))
          state.set(i, "Emp");
        for (auto i : range1(n_up))
          state.set(i, "Up");
        for (auto i : range1(n_dn))
          state.set(i+n_up, "Dn");
        auto psi0 = randomMPS(state);

        auto [energy, psi] = dmrg(H, psi0, sweeps, {"Silent", true});

        // Compute the entanglement entropies

        std::cout << "else if (qn == qn_tj({" << n_up << "," << n_dn << "}))\n{\n";
        std::cout << "svns = {";

        for (int b = 1; b < N; ++b) {
          psi.position(b);

          auto l = leftLinkIndex(psi, b);
          auto s = siteIndex(psi, b);
          auto [U, S, V] = svd(psi(b), {l, s});
          auto u = commonIndex(U, S);

          Real SvN = 0.;
          for (auto n : range1(dim(u))) {
            auto Sn = elt(S, n, n);
            auto p = sqr(Sn);
            if (p > 1E-12)
              SvN += -p * log(p);
          }
          // printfln("Across bond b=%d, SvN = %.10f", b, SvN);
          if (b < N - 1)
            std::cout << SvN << ",\n";
          else
            std::cout << SvN << "};\n";
        }
        std::cout << "e0 = " << energy << ";\n";
        std::cout << "}\n";
      }
    }

  return EXIT_SUCCESS;
}
