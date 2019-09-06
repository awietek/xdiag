#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <chrono>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "hubbardthermotpq.options.h"

namespace lila { LoggerMPI loggermpi; }

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::models::HubbardModelMPI;
  using namespace hydra::utils;
  using namespace hydra::operators;
  using namespace hydra::thermodynamics;
  using hydra::thermodynamics::detail::combine_thermodynamics;
  using lila::loggermpi;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  
  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string ensemble = "canonical";
  std::string temperaturefile;
  int np = -1;
  double mu = 0.;
  int seed = 1;
  double mintemperature = 0.01;
  double precision = 1e-12;
  int iters = 1000;
  int verbosity = 1;
  parse_cmdline(outfile, latticefile, couplingfile, ensemble, temperaturefile, 
		np, mu, seed, mintemperature, precision, iters, verbosity, 
		argc, argv);

  loggermpi.out(1, "Using {} MPI tasks\n", mpi_size);

  // Input checks
  check_if_file_exists(latticefile);
  check_if_file_exists(couplingfile);
  check_if_file_exists(temperaturefile);

  // Parse bondlist/couplings/temperatures from file
  auto bondlist = read_bondlist(latticefile);
  auto couplings = read_couplings(couplingfile);
  auto temperatures = lila::ReadVector<double>(temperaturefile);

  // Select qns to compute
  int n_sites = bondlist.n_sites();
  std::vector<hubbard_qn> qns;
  if (ensemble == "canonical")
    {
      if (np < 0) np = n_sites;
      for (int nup=0; nup<=np; ++nup)
	{
	  int ndown = np - nup;
	  qns.push_back({nup, ndown});
	}
    }
  else if (ensemble == "grandcanonical")
    {
      loggermpi.err("grandcanonical not implemented yet!");
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

  // HOld results for all temperatures and quantum numbers
  int n_temperatures = temperatures.size();
  int n_qns = (int)qns.size();
  std::vector<double> e0s_for_qn(n_qns);
  std::vector<std::vector<double>> partitions_for_qn(n_temperatures);
  std::vector<std::vector<double>> energies_for_qn(n_temperatures);
  std::vector<std::vector<double>> quad_moments_for_qn(n_temperatures);
  for (int i=0; i < (int)n_temperatures; ++i)
    {
      partitions_for_qn[i].resize(n_qns, 0);
      energies_for_qn[i].resize(n_qns, 0);
      quad_moments_for_qn[i].resize(n_qns, 0);
    }

  // Compute TPQ for all quantum numbers
  int qn_idx = 0;
  for (auto qn : qns)
    {
      loggermpi.out(1, "Running nup={}, ndown={}\n", 
		    qn.n_upspins, qn.n_downspins);

      // Create output files
      std::stringstream ss;
      ss << outfile 
	 << ".nup." << qn.n_upspins 
	 << ".ndown." << qn.n_downspins;
      std::ofstream of;
      of.open(ss.str(), std::ios::out);
      if(of.fail()) 
	{
	  loggermpi.err("Could not open outfile with filename [{}]!", ss.str());
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}

      ss << ".tmatrix";
      std::ofstream of_tmatrix;
      of_tmatrix.open(ss.str(), std::ios::out);
      if(of_tmatrix.fail()) 
	{
	  loggermpi.err("Could not open tmatrixfile with filename [{}]!", 
			ss.str());
	  MPI_Abort(MPI_COMM_WORLD, 1);
	}


      loggermpi.out(1, "Initializing (works only up to 32 sites)...\n");
      auto t1 = MPI_Wtime();
      auto model = HubbardModelMPI<double>(bondlist, couplings, qn);
      auto t2 = MPI_Wtime();
      loggermpi.out(1, "Done. time init: {} secs\n", t2 - t1);

      // Create normalized, normal distributed start vector
      loggermpi.out(1, "Creating start vector ...\n");
      t1 = MPI_Wtime();
      lila::normal_dist_t<double> dist(0., 1.);
      lila::normal_gen_t<double> gen(dist, seed);
      lila::VectorMPI<double> v0(model.local_dim());
      lila::Random(v0, gen);
      // LilaPrint(seed + qn_idx*19277);
      // LilaPrint(v0.vector_local());
      lila::Normalize(v0);
      t2 = MPI_Wtime();
      loggermpi.out(1, "Done. time start vec: {} secs\n", t2 - t1);

      // Run the TPQ algorithm
      loggermpi.out(1, "Running TPQ algorithm ...\n");
      t1 = MPI_Wtime();
      auto res = thermodynamics_tpq(model, v0, temperatures, 
				    mintemperature, precision, 
				    iters, lila::loggermpi);
      t2 = MPI_Wtime();
      
      loggermpi.out(1, "Done. #TPQ steps: {}, time TPQ: {} secs\n\n",
		    res.tmatrix.size(), t2 - t1);

      std::stringstream sss;
      sss << "Block, nup: " << qn.n_upspins 
	  << ", ndown: " << qn.n_downspins;
      write_thermodynamics_tpq_outfile(res, of, sss.str());
      write_thermodynamics_tpq_tmatrix(res, of_tmatrix, sss.str());

      for (int i=0; i<n_temperatures; ++i)
	{
	  partitions_for_qn[i][qn_idx] = res.partitions(i);
	  energies_for_qn[i][qn_idx] = res.energies(i);
	  quad_moments_for_qn[i][qn_idx] = res.quad_moments(i);
	}
      e0s_for_qn[qn_idx] = res.e0;

      ++qn_idx;
    }  // for (auto qn : qns)

  std::vector<double> tot_partitions, tot_energies, 
    tot_quad_moments, tot_specific_heats;
  double tot_e0 = *std::min_element(e0s_for_qn.begin(), e0s_for_qn.end());
  combine_thermodynamics(temperatures, e0s_for_qn, partitions_for_qn,
			 energies_for_qn, quad_moments_for_qn, 
			 tot_partitions, tot_energies, 
			 tot_quad_moments, tot_specific_heats);
  thermodynamics_tpq_result_t tot_results;

  tot_results.e0 = tot_e0;
  tot_results.temperatures = temperatures;
  tot_results.partitions = tot_partitions;
  tot_results.energies = tot_energies;
  tot_results.quad_moments = tot_quad_moments;
  
  std::ofstream of;
  of.open(outfile, std::ios::out);
  if(of.fail()) 
    {
      loggermpi.err("Could not open outfile with filename [{}]!", outfile);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  write_thermodynamics_tpq_outfile(tot_results, of);

  MPI_Finalize();
  return EXIT_SUCCESS;
}
