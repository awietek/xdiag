#include <cstdlib>
#include <vector>
#include <utility>

#include <lila/all.h>
#include <hydra/all.h>

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::models::HubbardModel;
  using hydra::utils::range;

  int n_sites = 8;
  hubbard_qn qn = {n_sites/2, n_sites/2}; 
  bool periodic = false;
  double t = 1.;
  double U = 4.;
  auto temperatures = {.1}; //lila::linspace<double>(0.1, 1, 10);

  // Set 1d chain lattice
  std::vector<std::pair<int, int>> hoppings;
  for (int i : range<>(n_sites-1))
    hoppings.push_back({i, i+1});
  if (periodic) hoppings.push_back({n_sites-1, 0});
      
  auto model = HubbardModel(n_sites, hoppings);
  auto hamilton = model.matrix(t, U, qn);
  // LilaPrint(hamilton);
  
  // lila::Matrix<double> hdag = Transpose(hamilton) ;
  // assert(lila::close<double>(hamilton, hdag));

  // LilaPrint(hamilton);
  auto eigs = lila::EigenvaluesH(hamilton); 
  double e0 = eigs(0);
 
  hydra::hilbertspaces::Hubbard<unsigned int> hs(n_sites, qn);
  printf("dim: %d, gs energy (ED): %f\n", hs.size(), e0);
  
  for (double t : temperatures)
    {
      double beta = 1./t;

      auto exp_eigs = eigs;
      lila::Map(exp_eigs, [t, e0](double& e) { e = std::exp(-(e-e0)/t); });
      double partition = lila::Sum(exp_eigs);
      double energy = lila::Dot(eigs, exp_eigs) / partition;
      printf("%f %f\n", t, energy);

      // Create Neel state
      using hydra::uint32;
      hydra::hilbertspaces::hubbard_state<uint32> neel_config = {0, 0};
      for (int i=0; i<n_sites; ++i)
	{
	  if (i%2==0) neel_config.downspins |= 1 << i;
	  else neel_config.upspins |= 1 << i;
	}
      std::cout << hydra::hilbertspaces::Print(n_sites, neel_config) << std::endl;
      hydra::indexing::IndexHubbard<hydra::indexing::IndexTable<hydra::hilbertspaces::Spinhalf<uint32>, uint32>> indexing(hs);
      int neel_idx = indexing.index(neel_config);
      lila::Vector<double> neel_state(hamilton.ncols());
      neel_state(neel_idx) = 1.;
      auto expH = hamilton;
      lila::ExpH(expH, -beta/2., 'U');
      auto tpq = Mult(expH, neel_state);
      tpq = tpq / lila::Norm(tpq);
      double neel_energy = lila::Dot(tpq, Mult(hamilton, tpq));

      auto asdf = Mult(Mult(expH, Mult(hamilton, expH)), neel_state);
      printf("neel energy: %f\n", neel_energy);


      // printf("T: %f, beta: %f, energy: %f, partition: %f\n", t, beta, energy/n_sites, partition);
    }
   
  return EXIT_SUCCESS;
}
