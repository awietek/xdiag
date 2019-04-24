#include <cassert>
#include <cmath>
#include <chrono>

#include <hydra/thermodynamics/thermodynamics_tpq.h>
#include <hydra/thermodynamics/thermodynamics_detail.h>
#include <hydra/utils/iochecks.h>
#include <hydra/models/hubbardmodel.h>

namespace hydra { namespace thermodynamics {

    template <class model_t>
    thermodynamics_tpq_result_t
    thermodynamics_tpq(const BondList& bondlist, 
		       const Couplings& couplings,
		       const std::vector<typename model_t::qn_t>& qns,
		       const std::vector<double>& temperatures, 
		       int seed, int iters, double precision, int neval)
    {
      // Allocate result data
      int n_temperatures = temperatures.size();
      int n_qns = qns.size();
      thermodynamics_tpq_result_t result;
      result.partitions_for_qn.resize(n_temperatures);
      result.energies_for_qn.resize(n_temperatures);
      result.quad_moments_for_qn.resize(n_temperatures);
      for (int i=0; i < (int)n_temperatures; ++i)
	{
	  result.partitions_for_qn[i].resize(n_qns, 0);
	  result.energies_for_qn[i].resize(n_qns, 0);
	  result.quad_moments_for_qn[i].resize(n_qns, 0);
	}
      for (int i=0; i < (int)n_temperatures; ++i)
	{
	  result.alphas.resize(n_qns);
	  result.betas.resize(n_qns);
	}

      int qn_idx = 0;
      for (auto qn : qns)
	{
	  // Compute T-Matrix using Lanczos algorithm
	  auto model = model_t(bondlist, couplings, qn); 
	  auto multiply = [&model](const lila::Vector<double>& v, 
				   lila::Vector<double>& w) {
	    using Clock = std::chrono::high_resolution_clock;
	    using secs = std::chrono::duration<float>;
	    static int iter=0;
	    auto t1 = Clock::now();
	    model.apply_hamiltonian(v, w);
	    auto t2 = Clock::now();
	    printf("dim: %ld, iter: %d, time MVM: %3.4f\n", model.dim(), iter, secs(t2-t1).count()); 
	  };
	  int64 dim = model.dim();
	  auto lzs = lila::Lanczos<double, decltype(multiply)>
	    (dim, seed + 1234321*qn_idx, iters, precision, neval, multiply);
	  lila::Vector<double> eigs = lzs.eigenvalues();
	  double e0 = eigs(0);
	  result.e0s.push_back(e0);
	  result.alphas[qn_idx] = lzs.alphas();
	  result.betas[qn_idx] = lzs.betas();
	  auto tmat = lzs.tmatrix();
	  auto tmatm = lila::Matrix<double>(tmat);
	  auto teig = Eigen(tmat);
	  auto teigs = teig.eigenvalues;
	  auto Q = teig.eigenvectors;
	  double shift = *std::min_element(teigs.begin(), teigs.end());

	  // Evaluate energy and quad moments at every temperature
	  int t_idx = 0;
	  for (double t : temperatures)
	    {
	      double beta = 1. / t;
	      auto diag = lila::Zeros<double>(tmat.size(), tmat.size());
	      for (int j = 0; j < tmat.size(); ++j)
		diag(j, j) = exp(-(beta / 2.) * (teigs(j) - shift));
	      auto Texp = Mult(Mult(Q, diag), Transpose(Q));
	      auto Texp0 = Texp.col(0);
		  
	      double partition = Dot(Texp0, Texp0) * dim;
	      double energy = Dot(Texp0, Mult(tmatm, Texp0)) * dim;
	      double quad = 
		Dot(Texp0, Mult(tmatm, Mult(tmatm, Texp0))) * dim;
		  
	      result.partitions_for_qn[t_idx][qn_idx] = partition;
	      result.energies_for_qn[t_idx][qn_idx] = energy;
	      result.quad_moments_for_qn[t_idx][qn_idx] = quad;
	      ++t_idx;
	    }

	  ++qn_idx;
	}  // for (auto qn : qns)
      detail::combine_thermodynamics
	(temperatures, result.e0s, result.partitions_for_qn,
	 result.energies_for_qn, result.quad_moments_for_qn,
	 result.partitions, result.energies, 
	 result.quad_moments, result.specific_heats);
      return result;
    }

    
    template <class qn_t>
    void write_thermodynamics_tpq
    (const std::vector<qn_t>& qns, 
     const std::vector<double>& temperatures,
     const thermodynamics_tpq_result_t& result,
     std::string outfile)
    {
      std::ofstream of;
      of.open(outfile, std::ios::out);
      utils::check_if_file_exists(outfile);
      of << "### method: exact\n";
      int qn_idx = 0;
      for (auto qn : qns)
	{
	  of << std::setprecision(20) << "## BLOCK: " << qn << "\n"	      
	     << "# gs energy: " << result.e0s[qn_idx] << "\n";
	  of << "# alphas betas\n";
	  for (int i=0; i<result.alphas[qn_idx].size(); ++i)
	    of << std::setprecision(20) 
	       << result.alphas[qn_idx](i) << " " 
	       << result.betas[qn_idx](i)<< "\n";
	  of << "# temperature partition energy quadmoment\n";      
	  int t_idx = 0;

	  for (double t : temperatures)
	    {
	      of << std::setprecision(20)  
		 << t << " " << result.partitions_for_qn[t_idx][qn_idx] << " " 
		 << result.energies_for_qn[t_idx][qn_idx] << " " 
		 << result.quad_moments_for_qn[t_idx][qn_idx] << "\n";
	      ++t_idx;
	    }
	  ++qn_idx;
	}
      of << "## COMBINATION\n";
      double total_e0 = *std::min_element(result.e0s.begin(), result.e0s.end());
      of << "# gs energy: " << total_e0 << "\n";
      of << "# temperature  partition  energy  quadmoments  specificheats\n";
      int t_idx = 0;
      for (double t : temperatures)
	{
	  of << std::setprecision(20)  
	     << t << " " 
	     << result.partitions[t_idx] << " " 
	     << result.energies[t_idx] << " " 
	     << result.quad_moments[t_idx] << " " 
	     << result.specific_heats[t_idx] << "\n";
	  ++t_idx;
	}
    }


    // HubbardModel
    using models::HubbardModel;

    template
    thermodynamics_tpq_result_t
    thermodynamics_tpq<HubbardModel>
    (const BondList& bondlist, const Couplings& couplings,
     const std::vector<HubbardModel::qn_t>& qns,
     const std::vector<double>& temperatures,
     int seed, int iters, double precision, int neval);

    template
    void write_thermodynamics_tpq<HubbardModel::qn_t>
    (const std::vector<HubbardModel::qn_t>& qns, 
     const std::vector<double>& temperatures,
     const thermodynamics_tpq_result_t& result,
     std::string outfile);

  }
}
