// Copyright 2019 Alexander Wietek - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef LILA_SPARSE_LOBPCG_H_
#define LILA_SPARSE_LOBPCG_H_

#include "../matrix.h"
#include "../vector.h"
#include "../common.h"
#include "gramschmidt.h"

namespace lila {
  
  template <class vector_t>
  struct LobpcgResults {
    using coeff_t = typename vector_t::coeff_type;
    Vector<real_t<coeff_t>> eigenvalues;
    std::vector<vector_t> eigenvectors;
    std::vector<Vector<real_t<coeff_t>>> eigenvalues_history;
  };
    
  template <class vector_t, class multiply_f>
  LobpcgResults<vector_t>
  Lobpcg(multiply_f A, std::vector<vector_t>& X, double tol, int maxiter)
  {
    using coeff_t = typename vector_t::coeff_type;
    
    // Init results and orthogonalize start vectors
    assert(has_full_rank(X));
    int m = (int)X.size();
    auto res_mus = Zeros<real_t<coeff_t>>(m);
    LobpcgResults<vector_t> res = {res_mus, gramschmidt(X), 
				   std::vector<Vector<real_t<coeff_t>>>()};

    auto& mus = res.eigenvalues;
    auto& Xb = res.eigenvectors;
    
    // Allocate vectors (TODO: avoid copy) 
    auto P = Xb;
    auto Pb = Xb;
    auto W = X;
    auto Wb = X;
    auto Xnew = Xb;
    auto Xbnew = Xb;
    auto Pnew = Xb;
    auto Pbnew = Xb;

    // Initialize X,Xb, and P
    for (int i = 0; i < m; ++i)
      {
	X[i] = Xb[i];
	A(X[i], Xb[i]);  // MVM	
	Zeros(P[i]);
	Zeros(Pb[i]);	
	Zeros(W[i]);
	Zeros(Wb[i]);
	Zeros(Xnew[i]);
	Zeros(Xbnew[i]);	
	Zeros(Pnew[i]);
	Zeros(Pbnew[i]);
      }    
  
    // Main iteration loop
    for (int iter=0; iter<maxiter; ++iter)
      {
	bool firstiter = (iter == 0);

	// Compute mus
	auto previous_mus = mus;
	for (int i = 0; i < m; ++i)
	  {
	    coeff_t mu = Dot(X[i], Xb[i]) / Dot(X[i], X[i]);
	    mus(i) = real(mu);
	    assert(close(imag(mu), (real_t<coeff_t>)0.));
	    W[i] = Xb[i] - mu * X[i];
	    A(W[i], Wb[i]);  // MVM
	  }
	res.eigenvalues_history.push_back(mus);

	// Check for convergence
	if (Norm(previous_mus - mus) < tol)
	  {
	    // printf("LOBPCG converged in %d steps\n", iter);
	    Xb = gramschmidt(Xb);  // TODO: orthonormalization throughout algorithm
	    break;  
	  }
	// Create gram matrices (skip P vectors on first iteration)
	int gramsize = (firstiter) ? 2*m : 3*m;	
	auto gramVAV = Zeros<coeff_t>(gramsize, gramsize);
	auto gramVV = Zeros<coeff_t>(gramsize, gramsize);
	std::vector<vector_t*> V;
	for (int i=0; i<m; ++i) V.push_back(&(W[i]));
	for (int i=0; i<m; ++i) V.push_back(&(X[i]));
	if (!firstiter) 
	  for (int i=0; i<m; ++i) V.push_back(&(P[i]));

	std::vector<vector_t*> Vb;
	for (int i=0; i<m; ++i) Vb.push_back(&(Wb[i]));
	for (int i=0; i<m; ++i) Vb.push_back(&(Xb[i]));
	if (!firstiter) 
	  for (int i=0; i<m; ++i) Vb.push_back(&(Pb[i]));

	for (int i=0; i<gramsize; ++i)
	  for (int j=i; j<gramsize; ++j)  // j=i: only top triangle 
	                                  //used by eigensolver
	    {
	      gramVAV(i, j) = Dot(*(V[i]), *(Vb[j]));
	      gramVV(i, j) = Dot(*(V[i]), *(V[j]));
	    }	    

	// Solve generalized eigenvalue problem
	auto evals = EigenGenSymDefInplace(gramVAV, gramVV, true, 'U', 1);
	auto& evecs = gramVAV;

	// Assemble new vectors
	for (int i=0; i<m; ++i)
	  {
	    auto evec = evecs.col(i);
	    Zeros(Xnew[i]);
	    Zeros(Xbnew[i]);
	    Zeros(Pnew[i]);
	    Zeros(Pbnew[i]);
	    for (int j=0; j<m; ++j)
	      {
		coeff_t alpha = evec(j);
		coeff_t beta = evec(j+m);
		coeff_t gamma = (firstiter) ? 0 : evec(j+2*m);
		Xnew[i] += alpha*W[j] + beta*X[j] + gamma*P[j];
		Xbnew[i] += alpha*Wb[j] + beta*Xb[j] + gamma*Pb[j];
		Pnew[i] += alpha*W[j] + gamma*P[j];
		Pbnew[i] += alpha*Wb[j] + gamma*Pb[j];		
	      }
	  }

	// Swap new and old vectors
	std::swap(X, Xnew);
	std::swap(Xb, Xbnew);
	std::swap(P, Pnew);
	std::swap(Pb, Pbnew);

      }  
    return res;
  }
  
}

#endif
