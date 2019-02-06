#include <lila/special.h>

#include <hydra/hilbertspaces/spinhalf.h>
#include <hydra/indexing/indextable.h>
#include <hydra/indexing/indexhubbard.h>
#include <hydra/utils/range.h>
#include <hydra/utils/bitops.h>

#include "hubbardmodel.h"

namespace hydra { namespace models {
    
    HubbardModel::HubbardModel(const int& n_sites,
			       const std::vector<std::pair<int, int>> neighbors)
      : n_sites_(n_sites),
	neighbors_(neighbors)
    {}

    std::vector<hilbertspaces::hubbard_qn> HubbardModel::quantumnumbers()
    {
      std::vector<hilbertspaces::hubbard_qn> qns;
      for (int n_upspins = 0; n_upspins <= n_sites_; ++n_upspins)
	for (int n_downspins = 0; n_downspins <= n_sites_; ++n_downspins)
	  qns.push_back({n_upspins, n_downspins});
      return qns;
    }

    lila::Matrix<double> HubbardModel::matrix
    (const double& t, const double& U, const hilbertspaces::hubbard_qn& qn) 
      const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbits;
      using utils::popcnt;
	
      Hubbard<uint32> hs(n_sites_, qn);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing(hs);
      int dim = indexing.size();
      lila::Matrix<double> hamilton(dim, dim);
      lila::Zeros(hamilton);

      // Apply Hubbard repulsion
      for (int idx : range<>(indexing.size()))
	{
	  auto state = indexing.state(idx);
	  hamilton(idx,idx) = U*(double)popcnt(state.upspins & state.downspins);
	}
      // Apply hopping terms
      for (auto pair : neighbors_)
	{
	  int s1 = std::min(pair.first, pair.second); 
	  int s2 = std::max(pair.first, pair.second);

	  uint32 flipmask = ((uint32)1 << s1) | ((uint32)1 << s2);

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
		  // printf("fermi %f\n", fermi);
		  auto new_state = state;
		  new_state.upspins = upspins ^ flipmask;
		  int new_idx = indexing.index(new_state);
		  hamilton(new_idx, idx) -= fermi*t;
		}

	      // downspins hopping
	      if (((downspins & flipmask) != 0) && 
		  ((downspins & flipmask) != flipmask))
		{
		  const double fermi = 
		    popcnt(gbits(downspins, s2-s1-1, s1+1)) % 2==0 ? 1. : -1.;
		  // printf("fermi %f\n", fermi);
		  auto new_state = state;
		  new_state.downspins = downspins ^ flipmask;
		  int new_idx = indexing.index(new_state);
		  hamilton(new_idx, idx) -= fermi*t;
		}
	    }
	}
      return hamilton;
    }


    lila::Vector<double> HubbardModel::apply_fermion
    (const lila::Vector<double>& state_before, hilbertspaces::hubbard_qn& qn, 
     std::string type, int site) 
      const
    {
      using hilbertspaces::Spinhalf;
      using hilbertspaces::Hubbard;
      using hilbertspaces::Print;
      using indexing::IndexTable;
      using indexing::IndexHubbard;
      using utils::range;
      using utils::gbit;
      using utils::gbits;
      using utils::popcnt;
	
      Hubbard<uint32> hs_before(n_sites_, qn);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing_before(hs_before);
      int dim_before = indexing_before.size();

      printf("type: %s\n", type.c_str());

      if (type == "cdagup") ++qn.n_upspins;
      else if (type == "cup") --qn.n_upspins;
      else if (type == "cdagdn") ++qn.n_downspins;
      else if (type == "cdn") --qn.n_downspins;
      else 
	{
	  std::cerr << "Error in apply_fermion: Invalid fermion type!" << std::endl;
	  exit(EXIT_FAILURE);
	}

      Hubbard<uint32> hs_after(n_sites_, qn);
      IndexHubbard<IndexTable<Spinhalf<uint32>, uint32>> indexing_after(hs_after);
      int dim_after = indexing_after.size();

      lila::Vector<double> state_after(dim_after);

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
		  // std::cout << "cdagup\n";
		  // std::cout << "before: " <<  Print(n_sites_, config) << "\n";
		  // std::cout << "after : " <<  Print(n_sites_, new_config) << "\n\n";
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(upspins, site-1, 0)) % 2==0 ? 1. : -1.;
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
		  // std::cout << "cup\n";
		  // std::cout << "before: " <<  Print(n_sites_, config) << "\n";
		  // std::cout << "after : " <<  Print(n_sites_, new_config) << "\n\n";
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(upspins, site-1, 0)) % 2==0 ? 1. : -1.;
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
		  // std::cout << "cdagdn\n";
		  // std::cout << "before: " <<  Print(n_sites_, config) << "\n";
		  // std::cout << "after : " <<  Print(n_sites_, new_config) << "\n\n";
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(downspins, site-1, 0)) % 2==0 ? 1. : -1.;
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
		  // std::cout << "cdn\n";
		  // std::cout << "before: " <<  Print(n_sites_, config) << "\n";
		  // std::cout << "after : " <<  Print(n_sites_, new_config) << "\n\n";
		  int new_idx = indexing_after.index(new_config);
		  double fermi = 
		    popcnt(gbits(downspins, site-1, 0)) % 2==0 ? 1. : -1.;
		  state_after(new_idx) += fermi * state_before(idx);
		}
	    }
	}
    
      return state_after;
    }
    
    
    
  }
};
