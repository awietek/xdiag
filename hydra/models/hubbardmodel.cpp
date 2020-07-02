#include <lila/all.h>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indextable.h>
#include <hydra/indexing/indexhubbard.h>
#include <hydra/utils/range.h>
#include <hydra/utils/bitops.h>

#include "hubbardmodel.h"
#include "hubbardmodeldetail.h"


namespace hydra { namespace models {
    
    template <class coeff_t>
    HubbardModel<coeff_t>::HubbardModel
    (BondList bondlist, Couplings couplings, hilbertspaces::hubbard_qn qn)
      : n_sites_(bondlist.n_sites()),
	qn_(qn)
    {
      hilbertspaces::Hubbard<uint32> hs(n_sites_, qn_);
      dim_ = hs.size();
      
      hubbardmodeldetail::set_hubbard_terms<coeff_t>
      (bondlist, couplings, hoppings_, hopping_amplitudes_,
       currents_, current_amplitudes_,
       interactions_, interaction_strengths_,
       onsites_, onsite_potentials_, 
       szszs_, szsz_amplitudes_,
       exchanges_, exchange_amplitudes_, U_);
    }

    template <class coeff_t>
    void HubbardModel<coeff_t>::set_qn(hilbertspaces::hubbard_qn qn) 
    { 
      qn_ = qn; 
      hilbertspaces::Hubbard<uint32> hs(n_sites_, qn_);
      dim_ = hs.size();
    }

    template <class coeff_t>
    lila::Matrix<coeff_t> HubbardModel<coeff_t>::matrix() const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;
	
      Hubbard<uint32> hs(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);
      int dim = indexing.size();
      lila::Matrix<coeff_t> hamilton(dim, dim);
      lila::Zeros(hamilton);

      // Apply Hubbard U term
      if (std::abs(U_) > 1e-14)
	{
	  for (int idx : range<>(indexing.size()))
	    {
	      auto state = indexing.state(idx);
	      hamilton(idx,idx) += U_*(double)popcnt(state.upspins & 
						     state.downspins);
	    }
	}

      // Apply Hubbard V term
      int interaction_idx=0;
      for (auto pair : interactions_)
	{
	  const int s1 = pair.first; 
	  const int s2 = pair.second;
	  const double V = interaction_strengths_[interaction_idx];
	  if (std::abs(V) > 1e-14)
	    {
	      for (int idx : range<>(indexing.size()))
		{
		  auto state = indexing.state(idx);
		  hamilton(idx,idx) += 
		    V * (double)((gbit(state.upspins, s1) + 
				  gbit(state.downspins, s1))*
				 (gbit(state.upspins, s2) + 
				  gbit(state.downspins, s2)));
		}
	    }
	  ++interaction_idx;
	}

      // Apply onsite chemical potential
      int onsite_idx=0;
      for (auto site : onsites_)
	{
	  const double mu = onsite_potentials_[onsite_idx];
	  if (std::abs(mu) > 1e-14)
	    {
	      for (int idx : range<>(indexing.size()))
		{
		  auto state = indexing.state(idx);
		  hamilton(idx,idx) -= 
		    mu * (double)((gbit(state.upspins, site) + 
				   gbit(state.downspins, site)));
		}
	    }
	  ++onsite_idx;
	}
      
      // Apply hopping terms
      int hopping_idx=0;
      for (auto pair : hoppings_)
	{
	  const int s1 = std::min(pair.first, pair.second); 
	  const int s2 = std::max(pair.first, pair.second);
	  const coeff_t t = hopping_amplitudes_[hopping_idx];
	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  const uint32 secondmask = (uint32)1 << s2;

	  if (std::abs(t) > 1e-14)
	    {
	      for (int idx : range<>(indexing.size()))
		  {
		  auto state = indexing.state(idx);
		  const uint32& upspins = state.upspins;
		  const uint32& downspins = state.downspins;

		  // upspins hopping
		  if (((upspins & flipmask) != 0) && 
		      ((upspins & flipmask) != flipmask))
		    {
		      const double fermi = 
			popcnt(gbits(upspins ^ downspins, s2-s1, s1)) & 1 ? 1. : -1.;
		      auto new_state = state;
		      new_state.upspins = upspins ^ flipmask;
		      int new_idx = indexing.index(new_state);
		      if (upspins & secondmask)
			hamilton(new_idx, idx) += fermi*t;
		      else 
			hamilton(new_idx, idx) -= lila::conj(t) * fermi;
       
		    }

		  // downspins hopping
		  if (((downspins & flipmask) != 0) && 
		      ((downspins & flipmask) != flipmask))
		    {
		      const double fermi = 
			popcnt(gbits(upspins ^ (downspins ^ secondmask), s2-s1, s1+1)) & 1 ? 1.:-1.;
		      auto new_state = state;
		      new_state.downspins = downspins ^ flipmask;
		      int new_idx = indexing.index(new_state);
		      if (downspins & secondmask)
		        hamilton(new_idx, idx) += fermi*t;
		      else 
			hamilton(new_idx, idx) -= lila::conj(t) * fermi;
		    }
		}
	    }
	  ++hopping_idx;
	}  // loop over hoppings

      // SzSz terms
      int szsz_idx = 0;
      for (auto pair : szszs_)
      {
        int s1 = pair.first;
        int s2 = pair.second;
        double jz = szsz_amplitudes_[szsz_idx]*0.25;

        if (std::abs(jz) > 1e-14)
        {
          for (int idx : range<>(indexing.size()))
          {
            auto state = indexing.state(idx);
            const uint32& upspins = state.upspins;
            const uint32& downspins = state.downspins;

            bool up1 = gbit(upspins, s1);
            bool up2 = gbit(upspins, s2);
            bool dn1 = gbit(downspins, s1);
            bool dn2 = gbit(downspins, s2);
//      			auto coeff = jz*(((double)gbit(upspins, s1) - (double)gbit(downspins, s1)) *
//			    ((double)gbit(upspins, s2) - (double)gbit(downspins, s2)));

            if (!(up1 && dn1) && !(dn2 && up2))
            {
              if ((up1 && up2) || (dn1 && dn2))
                hamilton(idx, idx) += jz;
              else if ((up1 && dn2) || (dn1 && up2))
                hamilton(idx, idx) += -jz;
          }
          }
        }
        ++szsz_idx;
      }

        // Exchange terms
        int exchange_idx=0;
        for (auto pair: exchanges_)
        {
          int s1 = std::min(pair.first, pair.second);
          int s2 = std::max(pair.first, pair.second);
          coeff_t jx = exchange_amplitudes_[exchange_idx]*0.5;
          const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
          if (std::abs(jx) > 1e-14)
            {
              for (int idx : range<>(indexing.size()))
          {
            auto state = indexing.state(idx);
            const uint32& upspins = state.upspins;
            const uint32& downspins = state.downspins;
          
      		  if ((popcnt(upspins & flipmask) == 1) &&
      			(popcnt(downspins & flipmask) == 1)&& popcnt((downspins & flipmask) & (upspins & flipmask)) == 0)
//if ((up1 && dn2 && !up2 && !dn1) || (dn1 && up2 && !dn2 && !up1))
              {
                auto new_state = state;
                new_state.upspins = upspins ^ flipmask;
                new_state.downspins = downspins ^ flipmask;

                if (popcnt(new_state.upspins) != popcnt(state.upspins))
            {
              using namespace hilbertspaces;
              std::cout << PrintSpinhalf(32, state.upspins)
                  << std::endl
                  << s1 << " " << s2
                  << std::endl
                  << popcnt(state.upspins) << std::endl
                  << PrintSpinhalf(32, new_state.upspins)
                  << std::endl
                  << popcnt(new_state.upspins) << std::endl
                  << std::endl;
            }
                assert(popcnt(new_state.upspins) == popcnt(state.upspins));
                assert(popcnt(new_state.downspins) == popcnt(state.downspins));

		            int flipped_idx = indexing.index(new_state);
                hamilton(flipped_idx, idx) += jx;
              }

          }  // loop over spin configurations
            }  // if (std::abs(jx) > 1e-14)
          ++exchange_idx;
        }  // for (auto pair: exchanges_)





      return hamilton;
    }

    template <class coeff_t>
    lila::Matrix<double> HubbardModel<coeff_t>::szMatrix(int siteIndex) const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;

      assert(0 <= siteIndex && siteIndex < n_sites_);
	
      Hubbard<uint32> hs(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);
      int dim = indexing.size();
      lila::Matrix<double> sz(dim, dim);
      lila::Zeros(sz);

      // Assemble sz matrix
      
        for (int idx : range<>(indexing.size()))
        {
          auto state = indexing.state(idx);
          const uint32& upspins = state.upspins;
          const uint32& downspins = state.downspins;
          bool ups = gbit(upspins, siteIndex);
          bool dns = gbit(downspins, siteIndex);

          sz(idx, idx) = ups - dns;
        }
      return sz;
      }  


    template <class coeff_t>
    lila::Matrix<double> HubbardModel<coeff_t>::sPlusMatrix(int siteIndex) const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;

      assert(0 <= siteIndex && siteIndex < n_sites_);
	
      Hubbard<uint32> hs(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);

      if (qn_.n_upspins < n_sites_ && qn_.n_downspins > 0) {

      // Initialize Hilbert space for target space
      auto target_qn = qn_;
      ++target_qn.n_upspins;
      --target_qn.n_downspins;
      Hubbard<uint32> target_hs(n_sites_, target_qn);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> target_indexing(target_hs);

      int base_dim = indexing.size();
      int target_dim= target_indexing.size();
      lila::Matrix<double> sPlus;
      sPlus.resize(target_dim, base_dim);
      lila::Zeros(sPlus);

      // Assemble sPlus matrix
    
      uint32 flipmask = ((uint32)1 << siteIndex);
      
        for (int idx : range<>(indexing.size()))
        {
          auto state = indexing.state(idx);
          const uint32& upspins = state.upspins;
          const uint32& downspins = state.downspins;
          bool ups = gbit(upspins, siteIndex);
          bool dns = gbit(downspins, siteIndex);

          if (dns && !ups) {
            auto flipped_state = state;
            flipped_state.upspins = state.upspins ^ flipmask;
            flipped_state.downspins = state.downspins ^ flipmask;
            int target_idx = target_indexing.index(flipped_state);

            sPlus(target_idx, idx) = 1;
          }
        }
      return sPlus;
      } else {
      lila::Matrix<double> sPlus(1, 1);
      lila::Zeros(sPlus);
      return sPlus;
      }
      }  


    template <class coeff_t>
    lila::Matrix<double> HubbardModel<coeff_t>::sMinusMatrix(int siteIndex) const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;

      assert(0 <= siteIndex && siteIndex < n_sites_);
	
      Hubbard<uint32> hs(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);

      if (qn_.n_downspins < n_sites_ && qn_.n_upspins > 0) {

      // Initialize Hilbert space for target space
      auto target_qn = qn_;
      --target_qn.n_upspins;
      ++target_qn.n_downspins;
      Hubbard<uint32> target_hs(n_sites_, target_qn);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> target_indexing(target_hs);

      int base_dim = indexing.size();
      int target_dim= target_indexing.size();
      lila::Matrix<double> sMinus;
      sMinus.resize(target_dim, base_dim);
      lila::Zeros(sMinus);

      // Assemble sMinus matrix
    
      uint32 flipmask = ((uint32)1 << siteIndex);
      
        for (int idx : range<>(indexing.size()))
        {
          auto state = indexing.state(idx);
          const uint32& upspins = state.upspins;
          const uint32& downspins = state.downspins;
          bool ups = gbit(upspins, siteIndex);
          bool dns = gbit(downspins, siteIndex);

          if (!dns && ups) {

            auto flipped_state = state;
            flipped_state.upspins = state.upspins ^ flipmask;
            flipped_state.downspins = state.downspins ^ flipmask;
            int target_idx = target_indexing.index(flipped_state);

            sMinus(target_idx, idx) = 1;
          }
        }
      return sMinus;
      } else {
      lila::Matrix<double> sMinus(1, 1);
      lila::Zeros(sMinus);
      return sMinus;
      }
      }  







    template <class coeff_t>
    void HubbardModel<coeff_t>::apply_hamiltonian
    (const lila::Vector<coeff_t>& in_vec, lila::Vector<coeff_t>& out_vec) const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;
	
      Hubbard<uint32> hs(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);
      int64 dim = indexing.size();
      assert(in_vec.size() == dim);
      assert(out_vec.size() == dim);
      lila::Zeros(out_vec);

      // Apply Hubbard U term
      if (std::abs(U_) > 1e-14)
	{
	  for (int idx : range<>(indexing.size()))
	    {
	      auto state = indexing.state(idx);
	      auto coeff =  U_*(double)popcnt(state.upspins & state.downspins);
	      out_vec(idx) += coeff * in_vec(idx);
	    }
	}

      // Apply Hubbard V term
      int interaction_idx=0;
      for (auto pair : interactions_)
	{
	  const int s1 = pair.first; 
	  const int s2 = pair.second;
	  const double V = interaction_strengths_[interaction_idx];
	  if (std::abs(V) > 1e-14)
	    {
	      for (int idx : range<>(indexing.size()))
		{
		  auto state = indexing.state(idx);
		  auto coeff = 	
		    V * (double)((gbit(state.upspins, s1) + 
				  gbit(state.downspins, s1))*
				 (gbit(state.upspins, s2) + 
				  gbit(state.downspins, s2)));
		  out_vec(idx) += coeff * in_vec(idx); 

		}
	    }
	  ++interaction_idx;
	}

      // Apply onsite chemical potential
      int onsite_idx=0;
      for (auto site : onsites_)
	{
	  const double mu = onsite_potentials_[onsite_idx];
	  if (std::abs(mu) > 1e-14)
	    {
	      for (int idx : range<>(indexing.size()))
		{
		  auto state = indexing.state(idx);
		  auto coeff = 
		    mu * (double)((gbit(state.upspins, site) + 
				   gbit(state.downspins, site)));
		  out_vec(idx) -= coeff * in_vec(idx); 
		}
	    }
	  ++onsite_idx;
	}
      
      // Apply hopping terms
      int hopping_idx=0;
      for (auto pair : hoppings_)
	{
	  const int s1 = std::min(pair.first, pair.second); 
	  const int s2 = std::max(pair.first, pair.second);
	  const coeff_t t = hopping_amplitudes_[hopping_idx];
	  const uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);
	  for (int idx : range<>(indexing.size()))
	    {
	      auto state = indexing.state(idx);
	      const uint32& upspins = state.upspins;
	      const uint32& downspins = state.downspins;
		
	      // upspins hopping
	      if (((upspins & flipmask) != 0) && 
		  ((upspins & flipmask) != flipmask))
		{
		  const double fermi = 
		    popcnt(gbits(upspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
		  auto new_state = state;
		  new_state.upspins = upspins ^ flipmask;
		  int new_idx = indexing.index(new_state);
		  out_vec(new_idx) -= fermi*t*in_vec(idx);
		}

	      // downspins hopping
	      if (((downspins & flipmask) != 0) && 
		  ((downspins & flipmask) != flipmask))
		{
		  const double fermi = 
		    popcnt(gbits(downspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
		  auto new_state = state;
		  new_state.downspins = downspins ^ flipmask;
		  int new_idx = indexing.index(new_state);
		  out_vec(new_idx) -= fermi*t*in_vec(idx);
		}
	    }
	  ++hopping_idx;
	}  // loop over hoppings
    }

    template <class coeff_t>
    hilbertspaces::hubbard_qn HubbardModel<coeff_t>::apply_fermion
    (const lila::Vector<coeff_t>& state_before, 
     lila::Vector<coeff_t>& state_after, 
     std::string type, int site) 
      const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;
	
      Hubbard<uint32> hs_before(n_sites_, qn_);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing_before(hs_before);

      auto qn_after = qn_;
      if (type == "cdagup") ++qn_after.n_upspins;
      else if (type == "cup") --qn_after.n_upspins;
      else if (type == "cdagdn") ++qn_after.n_downspins;
      else if (type == "cdn") --qn_after.n_downspins;
      else 
	{
	  std::cerr << "Error in apply_fermion: Invalid fermion type!" << std::endl;
	  exit(EXIT_FAILURE);
	}

      Hubbard<uint32> hs_after(n_sites_, qn_after);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing_after(hs_after);
      int64 dim_after = indexing_after.size();

      state_after.resize(dim_after);

      const uint32 sitemask = (1 << site);
      const uint32 antisitemask = ~(1 << site);
      if (type == "cdagup")
	{
	  for (int idx : range<>(indexing_before.size()))
	    {
	      auto config = indexing_before.state(idx);
	      const uint32& upspins = config.upspins;

	      // raise local site val if 0
	      if (gbit(upspins, site) == 0) 
		{
		  auto new_config = config;
		  new_config.upspins |= sitemask;
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(upspins, site, 0)) % 2==0 ? 1. : -1.;
		  state_after(new_idx) += fermi * state_before(idx);
		}
	    }
	} 
      else if (type == "cup")
	{
	  for (int idx : range<>(indexing_before.size()))
	    {
	      auto config = indexing_before.state(idx);
	      const uint32& upspins = config.upspins;
	
	      // lower local site val if 1
	      if (gbit(upspins, site) == 1) 
		{
		  auto new_config = config;
		  new_config.upspins &= antisitemask;
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(upspins, site, 0)) % 2==0 ? 1. : -1.;
		  state_after(new_idx) += fermi * state_before(idx);
		}
	    }
	}
      else if (type == "cdagdn")
	{
	  for (int idx : range<>(indexing_before.size()))
	    {
	      auto config = indexing_before.state(idx);
	      const uint32& downspins = config.downspins;
	
	      // raise local site val if 0
	      if (gbit(downspins, site) == 0) 
		{
		  auto new_config = config;
		  new_config.downspins |= sitemask;
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(downspins, site, 0)) % 2==0 ? 1. : -1.;
		  state_after(new_idx) += fermi * state_before(idx);
		}
	    }
	} 
      else if (type == "cdn")
	{
	  for (int idx : range<>(indexing_before.size()))
	    {
	      auto config = indexing_before.state(idx);
	      const uint32& downspins = config.downspins;
	
	      // lower local site val if 1
	      if (gbit(downspins, site) == 1) 
		{
		  auto new_config = config;
		  new_config.downspins &= antisitemask;
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(downspins, site, 0)) % 2==0 ? 1. : -1.;
		  state_after(new_idx) += fermi * state_before(idx);
		}
	    }
	}
      return qn_after;
    }
    
    template class HubbardModel<double>;
    template class HubbardModel<complex>;    
  }
}
