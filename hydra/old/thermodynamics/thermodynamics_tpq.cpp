#include <cassert>
#include <cmath>
#include <chrono>

#include <hydra/thermodynamics/thermodynamics_tpq.h>
#include <hydra/thermodynamics/thermodynamics_detail.h>
#include <hydra/utils/iochecks.h>
#include <hydra/models/hubbardmodel.h>
#include <hydra/models/hubbardmodelmpi.h>

#include <lila/all.h>

namespace hydra { namespace thermodynamics {

    template <class model_t, class vector_t, class logger_t>
    thermodynamics_tpq_result_t
    thermodynamics_tpq(model_t& model, vector_t& v0,
		       std::vector<double> temperatures,
		       double mintemperature, double precision,
		       int iters, logger_t& logger)
    {
      using coeff_type = typename vector_t::coeff_type;
      using real_type = lila::real_t<coeff_type>;
      using tmatrix_t = lila::Tmatrix<real_type>;
      
      // Define multiplication
      auto multiply = [&model, &logger](const vector_t& v, vector_t& w) 
	{
	  using Clock = std::chrono::high_resolution_clock;
	  using secs = std::chrono::duration<float>;
	  static int iter=0;
	  auto t1 = Clock::now();
	  model.apply_hamiltonian(v, w);
	  auto t2 = Clock::now();
	  logger.out(2, "dim: {}, iter: {}, time MVM: {}\n", 
		     model.dim(), iter, secs(t2-t1).count());
	};
      

      double max_invt = 1. / mintemperature;
      double max_invt_half = max_invt / 2.;
      double norm = lila::Norm(v0);

      // Set convergence criterion
      auto converged = 
	[precision, iters, max_invt_half, norm](const tmatrix_t& tmat, real_type beta) 
	{ 
	  return LanczosConvergedTimeEvolution(tmat, beta, max_invt_half, 
					       precision, iters, norm);
	};

      // Call the Lanczos algorithm
      auto lczs_result = Lanczos(multiply, v0, converged);

      auto tmat = lczs_result.tmatrix;
      auto tmatm = lila::Matrix<double>(tmat);
      auto teig = Eigen(tmat);
      auto teigs = teig.eigenvalues;
      auto Q = teig.eigenvectors;
      double e0 = *std::min_element(teigs.begin(), teigs.end());
      
      thermodynamics_tpq_result_t result;
      result.e0 = e0;
      result.tmatrix = tmat;
      result.eigenvalues = teigs;

      uint64 dim = model.dim();

      // Evaluate energy and quad moments at every temperature
      for (double t : temperatures)
	{
	  double beta = 1. / t;
	  auto diag = lila::Zeros<double>(tmat.size(), tmat.size());
	  for (int j = 0; j < tmat.size(); ++j)
	    diag(j, j) = exp(-(beta / 2.) * (teigs(j) - e0));
	  auto Texp = Mult(Mult(Q, diag), Transpose(Q));
	  auto Texp0 = Texp.col(0);
		  
	  double partition = Dot(Texp0, Texp0) * dim;
	  double energy = Dot(Texp0, Mult(tmatm, Texp0)) * dim;
	  double quad_moment = 
	    Dot(Texp0, Mult(tmatm, Mult(tmatm, Texp0))) * dim;
		  
	  result.temperatures.push_back(t);
	  result.partitions.push_back(partition);
	  result.energies.push_back(energy);
	  result.quad_moments.push_back(quad_moment);
	}

      return result;
    }


    void write_thermodynamics_tpq_outfile
    (const thermodynamics_tpq_result_t& res, std::ofstream& of, std::string comment)
    {
      if (comment != "")
	of << "## "<< comment << "\n";
      of << "# gs energy: " << std::setprecision(20) << res.e0 << "\n";
      of << "# temperature  partition  energy  quadmoments\n";
      int t_idx = 0;
      for (double t : res.temperatures)
	{
	  of << std::setprecision(20)  
	     << t << " " 
	     << res.partitions(t_idx) << " " 
	     << res.energies(t_idx) << " " 
	     << res.quad_moments(t_idx) << "\n"; 
	  ++t_idx;
	}
    }

    void write_thermodynamics_tpq_tmatrix
    (const thermodynamics_tpq_result_t& res, std::ofstream& of, std::string comment)
    {
      if (comment != "")
	of << "## "<< comment << "\n";
      of << "# tmatrix\n";
      of << "# alphas betas eigenvalues\n";
      auto alphas = res.tmatrix.diag();
      auto betas = res.tmatrix.offdiag();
      betas.push_back(0.);
      auto eigs = res.eigenvalues;
      assert(alphas.size() == betas.size());
      assert(alphas.size() == eigs.size());
      for (int i=0; i<alphas.size(); ++i)
	of << std::setprecision(20)  
	   << alphas(i) << " " << betas(i) << " " << eigs(i) << "\n"; 
    }



    template thermodynamics_tpq_result_t
    thermodynamics_tpq<hydra::models::HubbardModelMPI<double, uint32>, 
		       lila::VectorMPI<double>, lila::LoggerMPI>
    (hydra::models::HubbardModelMPI<double, uint32>& model, 
     lila::VectorMPI<double>& v0, std::vector<double> temperatures,
     double mintemperature, double precision, int iters, 
     lila::LoggerMPI& logger);


      
  }
}
