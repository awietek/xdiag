#include <cstdlib>

#include <lila/all.h>
#include <hydra/all.h>
#include <lime/all.h>

#include <stdio.h>

lila::Logger lg;

#include "hubbardfullonebodydensity.options.h"

template <class coeff_t>
void run_real_complex(std::string real_complex,
    hydra::all::BondList bondlist,
    hydra::all::Couplings couplings,
    hydra::all::hubbard_qn qn,
    int verbosity, int seed,
    std::string outfile)
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  FileH5 file;
  file = lime::FileH5(outfile, "w!");

  lg.out(1, "Creating {} Hubbard model for n_upspins={}, n_downspins={}...\n",
      real_complex, qn.n_upspins, qn.n_downspins);
  
  // Create Hubbard Hamiltonian 
  
  int n_sites = bondlist.n_sites();
  auto model = HubbardModel<coeff_t>(bondlist, couplings, qn);

  auto H = model.matrix();
  lg.out(1, "dim: {}\n", FormatWithCommas(model.dim())); 
  
  lg.out(1, "Diagonalizing...\n",
	     qn.n_upspins, qn.n_downspins);

  auto diagResults = EigenSym(H);

  file["Eigenvalues"] = diagResults.eigenvalues;
  file["GroundStateEnergy"] = diagResults.eigenvalues(0);
  file["Eigenvectors"] = diagResults.eigenvectors;
  file["Dimension"] = model.dim();

  lg.out(1, "E0: {}\n", diagResults.eigenvalues(0));



  auto groundstate = diagResults.eigenvectors.col(0);

  coeff_t norm = Dot(groundstate, groundstate);
  groundstate /= sqrt(norm);
  LilaPrint(groundstate);
  lg.out(1, "Done\n");
  lg.out(1, "Ground state energy: {}\n", diagResults.eigenvalues(0));


  // Generate one-body density matrix

  lila::Matrix<complex> singleParticleCorrelations;
  singleParticleCorrelations.resize(n_sites, n_sites);

  // First, generate a bondlist with all-to-all coupling
  
  std::vector<Bond> bonds;
  for (int i=0;i<n_sites;i++) {
    for (int j=i;j<n_sites;j++) {
      if (i != j)
      {
        std::vector<int> sites{i, j};
        std::string hopping_label = std::to_string(i) + "_" + std::to_string(j);
        bonds.push_back(Bond("HUBBARDHOP", hopping_label, sites));
      } else {
        std::vector<int> sites{i};
        std::string hopping_label = std::to_string(i);
        bonds.push_back(Bond("HUBBARDMU", hopping_label, sites));
      }
    }
  }
  BondList densityBondlist(bonds);
  
  // Measure diagonal elements first

  for (int i=0;i<n_sites;i++) {
    std::map<std::string, complex> coupling_map;
    std::string hopping_label = std::to_string(i);
    coupling_map[hopping_label] = -1;
    Couplings densityCouplings(coupling_map);
    auto hoppingModel = HubbardModel<coeff_t>(densityBondlist, densityCouplings, qn);
    auto hoppingH = hoppingModel.matrix();
    singleParticleCorrelations(i, i) = Dot(groundstate, Mult(hoppingH, groundstate));
    std::cout << i << std::endl;
    std::cout << singleParticleCorrelations(i,i) << std::endl;
  }

  // Then measure off-diagonal elements.
  // If hoppings are complex, it is necessary for the one-body density matrix
  // to measure <c_i^\dagger c_j> and <c_j^\dagger c_i> separately.
  // Note: the technical reason this exception is necessary is because complex
  // hoppings won't work if coeff_t is real, since the ground state will be real.
  // This would not be necessary if, in the future,
  // complex Hamiltonians can act on real ground states.
  if (real_complex == std::string("COMPLEX")) {
    for (int i=0;i<n_sites;i++) {
      for (int j=(i+1);j<n_sites;j++) {

        // First measure <c_i^\dagger c_j + c_j^\dagger c_j>
        
        std::map<std::string, complex> coupling_map;
        std::string hopping_label = std::to_string(i) + "_" + std::to_string(j);
        coupling_map[hopping_label] = 1;
        Couplings densityCouplings(coupling_map);
        auto hopping_symmetric = HubbardModel<coeff_t>(densityBondlist, densityCouplings, qn);
        auto hopping_symmetricH = hopping_symmetric.matrix();
        complex symmetric_expectation = Dot(groundstate, Mult(hopping_symmetricH, groundstate));

        // Then measure <c_i^\dagger c_j - c_j^\dagger c_i>
        coupling_map[hopping_label] = std::complex<double>(0, 1);
        densityCouplings = Couplings(coupling_map);
        auto hopping_antisymmetric = HubbardModel<coeff_t>(densityBondlist, densityCouplings, qn);
        auto hopping_antisymmetricH = hopping_antisymmetric.matrix();
        complex antisymmetric_expectation = complex(0, -1)*Dot(groundstate, Mult(hopping_antisymmetricH, groundstate));
        singleParticleCorrelations(i, j) = -0.5*(symmetric_expectation + antisymmetric_expectation);
        singleParticleCorrelations(j, i) = -0.5*(symmetric_expectation - antisymmetric_expectation);
      }
    }
  } else {
    for (int i=0;i<n_sites;i++) {
      for (int j=(i+1);j<n_sites;j++) {
        std::cout << i << " " << j << std::endl;
        std::map<std::string, complex> coupling_map;
        std::string hopping_label = std::to_string(i) + "_" + std::to_string(j);
        coupling_map[hopping_label] = 1;
        Couplings densityCouplings(coupling_map);
        auto hoppingModel = HubbardModel<coeff_t>(densityBondlist, densityCouplings, qn);
        auto hoppingH = hoppingModel.matrix();
        singleParticleCorrelations(i,j) = -0.5*Dot(groundstate, Mult(hoppingH, groundstate));
        std::cout << singleParticleCorrelations(i,j) << std::endl;
        singleParticleCorrelations(j,i) = singleParticleCorrelations(i,j);
      }
    }
  }

  file["Correlations"] = singleParticleCorrelations;
  file["Hoppings"] = model.single_particle_hopping();
	file.close();
}

int main(int argc, char* argv[])
{
  using namespace hydra::all;
  using namespace lila;
  using namespace lime;

  // Get input parameters
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  int nup = -1;
  int ndown = -1;
  int seed = 1;
  int verbosity = 1;
  parse_cmdline(outfile, latticefile, couplingfile, nup, ndown, 
		  verbosity, seed, argc, argv);
  lg.set_verbosity(verbosity);  
  check_if_files_exists({latticefile, couplingfile});


  // Parse bondlist, couplings, and corr_bondlist from file
  BondList bondlist = read_bondlist(latticefile);
  Couplings couplings = read_couplings(couplingfile);


  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};


  if (couplings.all_real())
    run_real_complex<double>(std::string("REAL"), bondlist,
        couplings, qn, verbosity, seed, outfile);
  else
    run_real_complex<lila::complex>(std::string("COMPLEX"), bondlist,
        couplings, qn, verbosity, seed, outfile);

  return EXIT_SUCCESS;
}
  
