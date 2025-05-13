// Copyright 2019 Alexander Wietek - All Rights Reserved.
// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

// partly derived from alps/ietl band lanczos method

#ifndef LILA_SPARSE_BANDLANCZOS_H_
#define LILA_SPARSE_BANDLANCZOS_H_

#include <utility>

#include "../matrix.h"
#include "../vector.h"
#include "../common.h"

namespace lila {

  namespace detail 
  {
    class indexer 
    {
    public:
      indexer(int pc)
	: location_(2*pc+1), 
	  deflated_(2*pc+1) 
      {
	pc_=pc;
	lastloc_=pc;
	lastvec_=pc;
	for (int i=0; i<(2*pc+1); ++i) 
	  {  
	    location_[i] = i;
	    deflated_[i] = false; 
	  }  
      }      

      int cnv(int j)
      {
	const auto it = std::find(location_.begin(),location_.end(),j);
	return std::distance(location_.begin(), it);
      }

      void next()
      {
	lastvec_++;
	do {
	  lastloc_++;
	  if (lastloc_ == (int)location_.size() )
	    lastloc_=0;
	} while (deflated_[lastloc_]);
	location_[lastloc_] = lastvec_;
      }

      void deflate(int old)
      {
	pc_--;
	const auto it = std::find(location_.begin(),location_.end(),old);
	int place = std::distance(location_.begin(), it);
	deflated_[place] = true;
      }  

    private:
      std::vector<int> location_;
      std::vector<bool> deflated_;
      int pc_;
      int lastloc_;
      int lastvec_;
    };
  }

  template <class coeff_t>
  struct lanczos_eigen_t 
  {
    Vector<real_t<coeff_t>> eigenvalues;
    std::vector<Vector<coeff_t>> eigenvectors;
    std::vector<int> multiplicity;
  };


  template <class coeff_t, class multiply_f, class vector_t = Vector<coeff_t>>
  class BandLanczos
  {

  public:
    BandLanczos(uint64 dimension, int random_seed, int max_iterations, 
		double precision, int num_eigenvalue, multiply_f multiply,
		int n_bands) 
      : dimension_(dimension),
	random_seed_(random_seed),
	max_iterations_(max_iterations),
	precision_(precision),
	num_eigenvalue_(num_eigenvalue),
	multiply_(multiply),
	n_bands_(n_bands),
	init_states_(nullptr)
    { }

    Matrix<coeff_t> tmatrix() const { return tmatrix_; }

    void set_init_states(std::vector<vector_t>& init_states)
    { 
      if (!has_full_rank(init_states))
	{
	  printf("Error in lila::BandLanczos: start states do not have full rank!\n");
	  exit(EXIT_FAILURE);
	}
      assert((int)init_states.size() == n_bands_);
      for (int k = 0; k < (int)init_states.size(); ++k)
	assert((uint64)init_states[k].size() == dimension_);
      init_states_ = &init_states; 
    }


    bool converged() const
    {
      auto tmat = tmatrix();
      int current_size = tmat.nrows();
      if ((current_size == 0) || (current_size == n_bands_)) return false;


      // Iteration count multible of number of bands
      if (current_size % n_bands_  == 0)
	{
	  auto prev_tmat = tmat;
	  int prev_size = current_size - n_bands_;
	  prev_tmat.resize(prev_size, prev_size);
	  auto eigs = EigenvaluesSym(tmat);
	  auto prev_eigs = EigenvaluesSym(prev_tmat);
	  double residue = std::abs((prev_eigs(num_eigenvalue_) - 
				     eigs(num_eigenvalue_)) / eigs(num_eigenvalue_));
	  return (residue < precision_);
	}
      else return false;
    }


    lanczos_eigen_t<coeff_t> eigenvalues(double deflation_tol = 1e-8)
    {
      using detail::indexer;

      int p = n_bands_;
      indexer index(p);
      int pc = p;
      std::vector<int> I;

      // Create and initialize band Lanczos vectors
      std::vector<vector_t> v(2*p + 1);
      lila::normal_dist_t<real_t<coeff_t>> ddist(0., 1.);
      lila::normal_gen_t<coeff_t> gen(ddist, random_seed_);
      for (int k=0; k < 2*p + 1; ++k)
	{
	  try 							       
	    {
	      // initialize first p vectors
	      if (k < p)
		{
		  if (!init_states_)
		    {
		      v[k].resize(dimension_);	      
		      Random(v[k], gen);
		    }
		  else 
		    {
		      assert((uint64)(*init_states_)[k].size() == dimension_);
		      std::swap(v[k], (*init_states_)[k]);
		    }
		}
	      else v[k].resize(dimension_);
	    }
	  catch(...)
	    {
	      std::cerr << "Error: Could not allocate Band Lanczos vector " << k 
			<< "!" << std::endl << std::flush;
	      exit(EXIT_FAILURE);
	    }
	}

      // Main Lanczos loop
      for (int j=0; ((j <= max_iterations_) && (pc > 0)); ++j) 
	{
	  auto norm = Norm(v[index.cnv(j)]);
	  tmatrix_.resize(j+1, j+1);
	  // Deflate if necessary
	  if (norm < deflation_tol)
	    {
	      printf("deflate\n");
	      if (j-pc >= 0)
		{
		  I.push_back(j-pc);
		  index.deflate(j-pc);
		}
	      pc--;
	      if (pc == 0) break;
	      else
		{
		  // set v_k <- v_{k+1}
		  for (int k=j; k < j+pc; ++k)
		    std::swap(v[index.cnv(k)], v[index.cnv(k+1)]); 

		  // return to compute |v_j|_2
		  --j;
		  continue;
		}
	    }

	  // Normalize vector
	  if (j-pc>=0)
	    tmatrix_(j,j-pc) = norm;
	  Scale((coeff_t)(1./norm), v[index.cnv(j)]);

	  // Orthogonalize against previous vectors
	  for (int k=j+1; k < j+pc; ++k) 
	    {
	      coeff_t alpha = Dot(v[index.cnv(j)], v[index.cnv(k)]);
	      if (k-pc>=0)
	        tmatrix_(j, k-pc) = alpha;
	      Add(v[index.cnv(j)], v[index.cnv(k)], -alpha);
	    }
	    
	  // Compute the matrix-vector product
	  multiply_(v[index.cnv(j)], v[index.cnv(j+pc)]);

	  // Orthogonalize next Lanczos vector
	  int k0 = std::max(j-pc, 0);
	  for (int k=k0; k<j; ++k) 
	    {               
	      tmatrix_(k,j) = lila::conj(tmatrix_(j,k));
	      Add(v[index.cnv(k)], v[index.cnv(j+pc)], -tmatrix_(k,j));
	    }

	  // Orthogonalize against deflated vectors
	  std::sort(I.begin(),I.end());
	  for (int i=0; i < (int)I.size()+1; ++i) 
	    {            
	      int k = (i==(int)I.size()) ? j : I[i];	     
	      coeff_t dot = Dot(v[index.cnv(k)], v[index.cnv(j+pc)]);
	      tmatrix_(k,j) = dot;
	      Add(v[index.cnv(k)], v[index.cnv(j+pc)], -dot);
	    }

	  for (int i=0; i<(int)I.size(); ++i) 
	    {               
	      int k=I[i];
	      tmatrix_(j,k) = lila::conj(tmatrix_(k,j));
	    }

	  index.next();

	  if (converged()) break;
	   
	
	}  // main Lanczos loop

      return eigenvalue_approx();
    }

    lanczos_eigen_t<coeff_t> eigenvalue_approx(int dim=-1, 
					       double ghost_tol = 1e-8, 
					       bool ghost_discarding=true)
    {
      assert(dim < tmatrix_.nrows());
      if (dim < 0) dim = tmatrix_.nrows();

      Vector<real_t<coeff_t>> eigenvals;
      std::vector<int> multiplicity;
            
      std::vector<int> mult(dim, 1);
      std::vector<int> ev_pos(dim, 0); // position of n_th eigenvector
      std::vector<bool> ghost(dim, false); 
 
      auto tmat_eigen = EigenSym(tmatrix_);
      auto tmat_eigs = tmat_eigen.eigenvalues;

      // Check for degeneracies and compute number of distinct eigenvalues
      int num_evs = 1;
      for (int i=0; i<dim-1; ++i) 
	{
	  // use ghost_tol as precision for degeneracy
	  if ( std::abs((tmat_eigs(i) - tmat_eigs(i+1)) ) > ghost_tol)
	    {
	      ev_pos[num_evs] = i+1;
	      mult[num_evs] = 1;
	      ghost[num_evs]= false;
	      num_evs++;
	    } 
	  else
	    mult[num_evs-1]++;
	}

      // Extract eigenvectors
      std::vector<Vector<coeff_t>> eigenvecs;
      for (int i=0; i<num_evs; ++i) 
	{
	  Vector<coeff_t> evec = tmatrix_.col(ev_pos[i]);
	  eigenvecs.push_back(evec);
	}

      // Ghost discarding, as described in Cullum, Willoughby (1981)
      int num_ghosts=0;
      if (ghost_discarding)
	{
	  // Build matrix T2 (first column and row deleted)
	  Matrix<coeff_t> T2(dim-1, dim-1);
	  for(int i=0; i<dim-1; ++i)
	    for(int j=0; j<dim-1; ++j)
	      T2(i, j) = tmatrix_(i+1, j+1);
	  Vector<real_t<coeff_t>> t2_eigs = EigenvaluesSym(T2);
      

	  // Check whether eigenvalues of T2 are close to T
	  // if so, discard as ghosts
	  for (int i=0; i<num_evs; ++i) 
	    {
	      if (mult[i]==1) 
		{
		  for (int j=0; j<dim-1; ++j) 
		    {
		      if (std::abs(t2_eigs(j)-tmat_eigs(ev_pos[i])) < ghost_tol)
			{
			  // printf("ghost found!\n");
			  ghost[i]=true;
			  eigenvecs.erase(eigenvecs.begin()+i-num_ghosts);
			  ++num_ghosts;
			}
		    }
		}
	    }
	}  // ghost_discarding

      // return only non-ghost eigenvalues
      for (int i=0; i<num_evs; ++i) 
	{
          if (!ghost[i]) 
	    {
	      eigenvals.push_back(tmat_eigs(ev_pos[i]));
	      multiplicity.push_back(mult[i]);
	    }
        }
      
      lanczos_eigen_t<coeff_t> res;
      res.eigenvalues = eigenvals;
      res.multiplicity = multiplicity;
      for (int i=0; i<(int)eigenvecs.size(); ++i) 
	res.eigenvectors.push_back(eigenvecs[i]);
      return res;
    }


  private:
    uint64 dimension_;
    int random_seed_;
    int max_iterations_;
    double precision_;
    int num_eigenvalue_;
    multiply_f multiply_;
    int n_bands_;
    Matrix<coeff_t> tmatrix_;

    std::vector<vector_t>* init_states_;

  };

      
}

#endif
