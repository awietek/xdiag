#include <algorithm>
#include <lila/all.h>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indextable.h>
#include <hydra/indexing/indexhubbard.h>
#include <hydra/utils/range.h>
#include <hydra/utils/bitops.h>

#include "tjmodel.h"
#include "hubbardmodeldetail.h"


namespace hydra { namespace models {
    
    template <class coeff_t>
    TJModel<coeff_t>::TJModel
    (BondList bondlist, Couplings couplings, hilbertspaces::hubbard_qn qn)
      : n_sites_(bondlist.n_sites()),
	qn_(qn)
    {
      using hydra::combinatorics::binomial;

      int nup = qn.n_upspins;
      int ndn = qn.n_downspins;
      assert(nup + ndn <= n_sites_);
      dim_ = binomial(n_sites_, nup) * binomial(n_sites_ - nup, ndn);

      // Currently unused operatuors
      std::vector<std::pair<int, int>> currents_;
      std::vector<coeff_t> current_amplitudes_;
      std::vector<std::pair<int, int>> interactions_;
      std::vector<double> interaction_strengths_;
      double U_;

      hubbardmodeldetail::set_hubbard_terms<coeff_t>
      (bondlist, couplings, hoppings_, hopping_amplitudes_,
       currents_, current_amplitudes_,
       interactions_, interaction_strengths_,
       onsites_, onsite_potentials_, 
       szszs_, szsz_amplitudes_,
       exchanges_, exchange_amplitudes_, U_);
    }

    template <class coeff_t>
    lila::Matrix<coeff_t> TJModel<coeff_t>::matrix(bool ninj_term) const
    {
      using state_t = uint32;
      using hydra::combinatorics::up_hole_to_down;
      using hydra::hilbertspaces::Spinhalf;
      using namespace hydra::utils;


      int nup = qn_.n_upspins;
      int ndn = qn_.n_downspins;
      
      // Try allocating upspin / downspin vectors
      std::vector<state_t> upspins;
      std::vector<state_t> dnspins;
      try
	{
	  upspins.resize(dim_);
	  dnspins.resize(dim_);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate upspin/downspin " 
		    << "vectors for TJModel!" << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // Fill vectors holding upspin/downspin configurations
      int64 idx=0;
      auto hs_upspins = Spinhalf<state_t>(n_sites_, nup);
      auto hs_holes_in_ups = Spinhalf<state_t>(n_sites_ - nup, ndn);
      for (state_t ups : hs_upspins)
	for (state_t holes : hs_holes_in_ups)
	  {
	    state_t dns = up_hole_to_down(ups, holes);
	    upspins[idx] = ups;
	    dnspins[idx] = dns;
	    ++idx;
	  }
      assert(idx == dim_);


      // Define lambda function to get indices
      auto index_of_up_dn = 
	[&upspins, &dnspins](state_t const& ups, state_t const& dns)
	{
	  // Binary search the new indices
	  auto up_bounds = std::equal_range(upspins.begin(), 
	  				    upspins.end(), 
	  				    ups);
	  auto up_begin = std::distance(upspins.begin(),
					up_bounds.first); 
	  auto up_end = std::distance(upspins.begin(),
				      up_bounds.second); 
	  auto it = 
	    std::lower_bound(dnspins.begin() + up_begin, 
			     dnspins.begin() + up_end, dns);
	  return std::distance(dnspins.begin(), it);
	  
	  // auto idxb = std::distance(dnspins.begin(), it);
	  // // linear search
	  // int64 idx = 0;
	  // for (idx=0; idx<(int64)upspins.size(); ++idx)
	  //   {
	  //     if ((ups == upspins[idx]) && (dns == dnspins[idx]))
	  // 	break;
	  //   }
	  // printf("linear: %d, binary: %d\n", idx, idxb);
	  // if (idx != idxb)
	  //   {
	  //     std::cout << "up_begin " << up_begin << ", up_end " << up_end << std::endl;
	  //   }
	  // assert(idx == idxb);
	  // return idx;
	};

      // Try allocating the matrix
      lila::Matrix<coeff_t> H;
      try 
	{
	  H.resize(dim_, dim_);
	  Zeros(H);
	}
      catch(...)
	{
	  std::cerr << "Error: Could not allocate matrix for TJModel!" 
		    << std::endl << std::flush;
	  exit(EXIT_FAILURE);
	}

      // SzSz terms
      int szsz_idx = 0;
      for (auto pair : szszs_)
	{
	  int s1 = pair.first;
	  int s2 = pair.second;
	  double jz = szsz_amplitudes_[szsz_idx]*0.25;

	  if (std::abs(jz) > 1e-14)
	    {
	      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = upspins[idx];
		  state_t dns = dnspins[idx];

		  bool up1 = gbit(ups, s1);
		  bool up2 = gbit(ups, s2);
		  bool dn1 = gbit(dns, s1);
		  bool dn2 = gbit(dns, s2);

		  if (ninj_term)
		    {
		      if ((up1 && up2) || (dn1 && dn2))
			H(idx, idx) += 0;
		      else if ((up1 && dn2) || (dn1 && up2))
			H(idx, idx) += -2*jz;
		    }
		  else
		    {
		      if ((up1 && up2) || (dn1 && dn2))
			H(idx, idx) += jz;
		      else if ((up1 && dn2) || (dn1 && up2))
			H(idx, idx) += -jz;
		    }
		} 

	    }  // if (std::abs(jz) > 1e-14)
	  ++szsz_idx;
	}  // for (auto pair : szszs_)

      // Exchange terms
      int exchange_idx=0;
      for (auto pair: exchanges_)
	{
	  int s1 = std::min(pair.first, pair.second);
	  int s2 = std::max(pair.first, pair.second);
	  coeff_t jx = exchange_amplitudes_[exchange_idx]*0.5;
	  state_t flipmask = ((state_t)1 << s1) | ((state_t)1 << s2);
	
	  if (std::abs(jx) > 1e-14)
	    {
	      for (int64 idx=0; idx<dim_; ++idx)
		{
		  state_t ups = upspins[idx];
		  state_t dns = dnspins[idx];
		  bool up1 = gbit(ups, s1);
		  bool up2 = gbit(ups, s2);
		  bool dn1 = gbit(dns, s1);
		  bool dn2 = gbit(dns, s2);
		
		  if ((up1 && dn2) || (dn1 && up2))
		    {
		      state_t flipped_ups = ups ^ flipmask;
		      state_t flipped_dns = dns ^ flipmask;
		      int64 flipped_idx = index_of_up_dn(flipped_ups, flipped_dns);
		      double fermi_up = 
			popcnt(gbits(ups, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
		      double fermi_dn = 
			popcnt(gbits(dns, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;

		      H(flipped_idx, idx) -= jx * fermi_up * fermi_dn;
		    }

		}  // loop over spin configurations
	    }  // if (std::abs(jz) > 1e-14)
	  ++exchange_idx;
	}  // for (auto pair: exchanges_)


      // Hoppings
      int hopping_idx=0;
      for (auto pair : hoppings_)
      	{
      	  int s1 = std::min(pair.first, pair.second); 
      	  int s2 = std::max(pair.first, pair.second);
      	  coeff_t t = hopping_amplitudes_[hopping_idx];
      	  uint32 flipmask = ((state_t)1 << s1) | ((state_t)1 << s2);
      	  if (std::abs(t) > 1e-14)
      	    {
      	      for (int64 idx=0; idx<dim_; ++idx)
      		{
      		  state_t ups = upspins[idx];
      		  state_t dns = dnspins[idx];
	      
      		  // upspin hopping
      		  if ((dns & flipmask) ==0)
      		    {
      		      if (((ups & flipmask) != 0) && ((ups & flipmask) != flipmask))
      			{
      			  state_t flipped_ups = ups ^ flipmask;
      			  double fermi_up = 
      			    popcnt(gbits(ups, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
      			  int64 flipped_idx = index_of_up_dn(flipped_ups, dns);
			  // printf("idx: %d, ups: %d, s1: %d, s2: %d, fups: %d, fidx: %d\n",
			  // 	 idx, ups, s1,s2, flipped_ups, flipped_idx);

      			  H(flipped_idx, idx) -= t * fermi_up;
      			}
      		    }

      		  // dnspin hopping
      		  if ((ups & flipmask) ==0)
      		    {
      		      if (((dns & flipmask) != 0) && ((dns & flipmask) != flipmask))
      			{
      			  state_t flipped_dns = dns ^ flipmask;
      			  double fermi_dn = 
      			    popcnt(gbits(dns, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
      			  int64 flipped_idx = index_of_up_dn(ups, flipped_dns);
      			  H(flipped_idx, idx) -= t * fermi_dn;
			  // printf("idx: %d, dns: %d, s1: %d, s2: %d, fdns: %d, fidx: %d\n",
			  // 	 idx, dns, s1,s2, flipped_dns, flipped_idx);
      			}
      		    }
      		}  // loop over spin configurations
      	    }  // if (std::abs(t) > 1e-14)
      	  ++hopping_idx;
      	}  // for (auto pair : hoppings_)
    
      return H;

    }

    template class TJModel<double>;
    template class TJModel<complex>;

  }
}
