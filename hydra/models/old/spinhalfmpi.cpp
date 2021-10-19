#include "spinhalfmpi.h"

#include <algorithm>
#include <mpi.h>

SpinhalfMPI::SpinhalfMPI()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size_);
}

inline int get_mpi_rank_spinhalf(uint64 const& spins, int const& mpi_size)
{ 
  uint64 x = spins;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x % mpi_size_;
}

void build_spinhalf_tables
(int mpi_rank, int mpi_size,
 std::pair<int,int> n_sites_n_upspins, int n_prefix,
 std::vector<LinTable>& lintables;
 std::map<std::pair<int,int>, std::vector<uint64>>& my_prefixes_,
 std::map<std::pair<int,int>, std::vector<uint64>>& my_prefix_size_, 
 std::map<std::pair<int,int>, std::vector<uint64>>& my_prefix_offset_)
{
  my_prefixes_[n_sites_n_upspins] = std::vector<uint64>();
  my_prefix_size_[n_sites_n_upspins] = std::vector<uint64>();
  my_prefix_offset_[n_sites_n_upspins] = std::vector<uint64>();
  auto& my_prefixes = my_prefixes_[n_sites_n_upspins];
  auto& my_prefix_size = my_prefix_size_[n_sites_n_upspins];
  auto& my_prefix_offset = my_prefix_offset_[n_sites_n_upspins];
  
  uint64 offset = 0;
  for (uint64 prefix=(uint64)0; prefix < (uint64)1<<n_prefix; ++prefix)
    {
      int n_upspins_prefix = popcnt(prefix);
      if (n_upspins_prefix > n_upspins) continue;

      if (get_mpi_rank_spinhalf(prefix, mpi_size) == mpi_rank)
	{
	  int n_upspins_postfix = n_upspins - n_upspins_prefix;
	  my_prefixes.push_back(prefix);
	  uint64 size = lintables[n_upspins_postfix].size();
	  my_prefix_size.push_back(size);
	  my_prefix_offset.push_back(offset);
	  offset += size;
	}
    }

}

uint64 SpinhalfMPI::dim(int n_sites, quantumnumber qn) const
{
  int n_upspins = qn_to_n_upspins(qn);
  return combinatorics::binomial(n_sites, n_upspins);
}

bool SpinhalfMPI::block_defined(int n_sites, quantumnumber qn)
{
  int n_sites = bondlist.n_sites();
  int n_upspins = qn_to_n_upspins(qn);
  std::pair<int,int> n_sites_n_upspins({n_sites, n_upspins});
  return (std::find(n_sites_n_upspins_.begin(), n_sites_n_upspins_.end(),
		    n_sites_n_upspins) == n_sites_n_upspins_.end())
}

void SpinhalfMPI::create_block(int n_sites, quantumnumber qn)
{
  // Create all Lin Tables
  if (std::find(n_sites_.begin(), n_sites_.end(), n_sites) == n_sites_.end())
    {
      lintables_[n_sites] = std::vector<LinTable>();
      for (int n_upspins=0; n_upspins <= n_sites; ++n_upspins)
	{
	  SpinhalfBasis basis(n_sites, n_upspins);
	  lintables_[n_sites].push_back(LinTable(basis));
	}
    }

  // Create prefix/postfix tables
  int n_upspins = qn_to_n_upspins(qn);
  int n_prefix = n_sites / 2;
  int n_postfix = n_sites - n_prefix;
  std::pair<int,int> n_sites_n_upspins({n_sites, n_upspins});
  if (std::find(n_sites_n_upspins_.begin(), n_sites_n_upspins_.end(),
		n_sites_n_upspins) == n_sites_n_upspins_.end())
    {
      n_sites_n_upspins_.push_back(n_sites_n_upspins);
      n_prefix_[n_sites_n_upspins] = n_prefix;
      n_postfix_[n_sites_n_upspins] = n_postfix;
      build_spinhalf_tables(mpi_rank_, mpi_size_, n_sites_n_upspins,
			    n_prefix, lintables,
			    my_prefixes, my_prefix_size, my_prefix_offset);
      build_spinhalf_tables(mpi_rank_, mpi_size_, n_sites_n_upspins,
			    n_postfix, lintables,
			    my_postfixes, my_postfix_size,my_postfix_offset);
      dim_local_[n_sites_n_upspins] =
	my_prefix_offset.back() + my_prefix_size.back();
    }
}

quantumnumber SpinhalfMPI::multiply(quantumnumber qn,
				    BondList const& bondlist,
				    Couplings const& couplings,
				    lila::VectorMPI<coeff_t> const& x,
				    lila::VectorMPI<coeff_t> & y,
				    coeff_t alpha = (coeff_t)1.0,
				    coeff_t beta  = (coeff_t)0.0)
{
  int n_sites = bondlist.n_sites();
  if (!block_defined(n_sites, qn))
    create_block(n_sites, qn);
  int n_upspins = qn_to_n_upspins(n_sites, qn);

  lila::Scale(beta, y);
  
  // Ising bonds
  auto ising_bonds =
    bondlist.bonds_of_type("Ising") +
    bondlist.bonds_of_type("Heisenberg") +
    bondlist.bonds_of_type("ISING") +
    bondlist.bonds_of_type("HEISENBERG");
  std::vector<std::tuple<int,int,double>> ising_s1_s2_J;
  for (auto const& bond : ising_bonds)
    {
      assert(bond.n_sites() == 2);
      int s1 = bond.sites(0);
      int s2 = bond.sites(1);
      assert(0 <= s1 < n_sites);
      assert(0 <= s2 < n_sites);
      auto coupling = bond.coupling();
      if (couplings.defined(coupling))
	{
	  double J = couplings[coupling];
	  if (std::abs(J)>1e-12) ising_s1_s2_J.push_back({s1, s2, J});
	}
    }
  multiplyIsing(n_sites, n_upspins, ising_s1_s2_J, x, y, alpha);
  
  // Exchange bonds
  auto exchange_bonds =
    bondlist.bonds_of_type("Exchange") +
    bondlist.bonds_of_type("Heisenberg") +
    bondlist.bonds_of_type("EXCHANGE") +
    bondlist.bonds_of_type("HEISENBERG");
  std::vector<std::tuple<int,int,double>> exchange_s1_s2_J;
  for (auto const& bond : ising_bonds)
    {
      assert(bond.n_sites() == 2);
      int s1 = bond.sites(0);
      int s2 = bond.sites(1);
      assert(0 <= s1 < n_sites);
      assert(0 <= s2 < n_sites);
      auto coupling = bond.coupling();
      if (couplings.defined(coupling))
	{
	  double J = couplings[coupling];
	  if (std::abs(J)>1e-12) exchange_s1_s2_J.push_back({s1, s2, J});
	}
    }
  multiplyExchange(n_sites, n_upspins, exchange_s1_s2_J, x, y, alpha);

}

void SpinhalfMPI::multiplyIsing(int n_sites, int n_upspins,
				std::vector<std::tuple<int,int,double>>
				const& s1_s2_Js,
				lila::VectorMPI<coeff_t> const& x,
				lila::VectorMPI<coeff_t> & y,
				coeff_t alpha = (coeff_t)1.0) const
{
  std::pair<int,int> n_sites_n_upspins({n_sites, n_upspins});

  uint64 dim_local = dim_local_[n_sites_n_upspins];
  int n_prefix = n_prefix_.at(n_sites_n_upspins);
  int n_postfix = n_postfix_.at(n_sites_n_upspins);
  auto& my_prefixes = my_prefixes_.at(n_sites_n_upspins);
  auto& my_prefix_size = my_prefix_size_.at(n_sites_n_upspins);
  auto& my_prefix_offset = my_prefix_offset_.at(n_sites_n_upspins);
  auto& my_postfixes = my_postfixes_.at(n_sites_n_upspins);
  auto& my_postfix_size = my_postfix_size_.at(n_sites_n_upspins);
  auto& my_postfix_offset = my_postfix_offset_.at(n_sites_n_upspins);
  
  for (auto s1_s2_J : s1_s2_Js)
    {
      int s1t, s2t;
      double J;
      std::tie(s1t, s2t, J) = s1_s2_J;
      int s1 = std::max(s1t, s2t); 
      int s2 = std::min(s1t, s2t); 
      
      uint64 s1_mask = (uint64)1 << s1;
      uint64 s2_mask = (uint64)1 << s2;
      for (uint64 prefix_idx=0; prefix_idx<(uint64)my_prefixes.size(),
	     ++prefix_idx)
	{
	  uint64 prefix = my_prefixes_[prefix_idx];
	  uint64 prefix_shift = prefix << n_prefix;
	  uint64 prefix_size = my_prefix_size[prefix_idx];
	  uint64 prefix_offset = my_prefix_size[prefix_offset];

	  int prefix_n_upspins = popcnt(prefix);
	  int postfix_n_upspins = n_upspins - prefix_n_upspins;
	  assert(postfix_n_upspins >= 0);

	  LinTable& lintable = lintables_[n_sites][postfix_n_upspins];
	  assert(lintable.size() == prefix_size);

	  uint64 begin = prefix_offset;
	  uint64 end = prefix_offset + prefix_size;
	  
	  // Case: s1 and s2 are both prefixes
	  if ((s1 >= n_postfix) && (s2 >= n_postfix))
	    {
	      coeff_t Jmult =  alpha * J / 4 ;
	      if (bool(prefix_shift | s1_mask) !=
		  bool(prefix_shift | s2_mask))
		double Jmult *= -1.;
	      for (uint64 idx=begin; idx<end; ++idx)
		y[idx] += x[idx] * Jmult;
	    }
	  
	  // Case: s1 is prefix, s2 is postfix
	  else if (s1 >= n_postfix)
	    {
	      bool s1_set = bool(prefix_shift | s1_mask);
	      coeff_t Jmult = alpha * J / 4;
	      for (uint64 idx=begin; idx<end; ++idx)
		{
		  uint64 postfix = lintable.state(idx);
		  y[idx] += x[idx] *
		    ((s1_set == bool(postfix | s2_mask)) ? Jmult : -Jmult);
		}
	    }

	  // Case: both s1 and s2 are postfixes
	  else
	    {
	      coeff_t Jmult = alpha * J / 4;
	      for (uint64 idx=begin; idx<end; ++idx)
		{
		  uint64 postfix = lintable.state(idx);
		  y[idx] += x[idx] *
		    ((bool(postfix | s1_mask) ==
		      bool(postfix | s2_mask)) ? Jmult : -Jmult);
		}
	    }
	  
	}  // loop over prefix indices   
    }  // loop over s1, s2, J
}


void SpinhalfMPI::multiplyExchange(int n_sites, int n_upspins,
				   std::vector<std::tuple<int,int,double>>
				   const& s1_s2_Js,
				   lila::VectorMPI<coeff_t> const& x,
				   lila::VectorMPI<coeff_t> & y,
				   coeff_t alpha = (coeff_t)1.0) const
{
  std::pair<int,int> n_sites_n_upspins({n_sites, n_upspins});

  uint64 dim_local = dim_local_[n_sites_n_upspins];
  int n_prefix = n_prefix_.at(n_sites_n_upspins);
  int n_postfix = n_postfix_.at(n_sites_n_upspins);
  auto& my_prefixes = my_prefixes_.at(n_sites_n_upspins);
  auto& my_prefix_size = my_prefix_size_.at(n_sites_n_upspins);
  auto& my_prefix_offset = my_prefix_offset_.at(n_sites_n_upspins);
  auto& my_postfixes = my_postfixes_.at(n_sites_n_upspins);
  auto& my_postfix_size = my_postfix_size_.at(n_sites_n_upspins);
  auto& my_postfix_offset = my_postfix_offset_.at(n_sites_n_upspins);

  // Find out which bonds only act on prefix/postfix sites or are mixed
  std::vector<std::tuple<int,int,double>> s1_s2_Js_postfix_only;
  std::vector<std::tuple<int,int,double>> s1_s2_Js_mixed;
  std::vector<std::tuple<int,int,double>> s1_s2_Js_prefix_only;
  for (auto s1_s2_J : s1_s2_Js)
    {
      int s1, s2;
      std::tie(s1, s2, std::ignore) = s1_s2_J;
      if ((s1 >= n_postfix) && (s2 >= n_postfix))
	s1_s2_Js_prefix_only.append(s1_s2_J);
      else if ((s1 < n_postfix) && (s2 < n_postfix))
	s1_s2_Js_postfix_only.append(s1_s2_J);
      else 
	s1_s2_Js_mixed.append(s1_s2_J);
    }
  
  // apply bonds which are only on the postfixes
  for (auto s1_s2_J : s1_s2_Js_postfix_only)
    {
      int s1t, s2t;
      double J;
      std::tie(s1t, s2t, J) = s1_s2_J;
      int s1 = std::max(s1t, s2t); 
      int s2 = std::min(s1t, s2t);
      
      uint64 flipmask = ((uint64)1 << s1) | ((uint64)1 << s1);
      coeff_t Jmult = alpha * J / 2;
      
      for (uint64 prefix_idx=0; prefix_idx<(uint64)my_prefixes.size(),
	     ++prefix_idx)
	{
	  uint64 prefix_size = my_prefix_size[prefix_idx];
	  uint64 prefix_offset = my_prefix_size[prefix_offset];

	  int prefix_n_upspins = popcnt(prefix);
	  int postfix_n_upspins = n_upspins - prefix_n_upspins;
	  assert(postfix_n_upspins >= 0);

	  LinTable& lintable = lintables_[n_sites][postfix_n_upspins];
	  assert(lintable.size() == prefix_size);

	  uint64 begin = prefix_offset;
	  uint64 end = prefix_offset + prefix_size;
	  for (uint64 idx=begin; idx<end; ++idx)
	    {
	      uint64 postfix = lintable.state(idx);
	      if (popcnt(postfix & flipmask)==1)
		{
		  uint64 postfix_flipped = postfix ^ flipmask;
		  uint64 new_idx = begin + lintable.index(postfix_flipped);
		  y[new_idx] += Jmult * x[idx]
		}
	    }
	  
	}  // loop over prefix indices   
    }  // loop over s1, s2, J for postfixes

  
  // apply bonds which act on both prefix and postfix
  for (auto s1_s2_J : s1_s2_Js_mixed)
    {
      int s1t, s2t;
      double J;
      std::tie(s1t, s2t, J) = s1_s2_J;
      int s1 = std::max(s1t, s2t); 
      int s2 = std::min(s1t, s2t);
      uint64 s1_mask = (uint64)1 << s1; 
      uint64 s2_mask = (uint64)1 << s2;
      
      coeff_t Jmult = alpha * J / 2;
      
      for (uint64 prefix_idx=0; prefix_idx<(uint64)my_prefixes.size(),
	     ++prefix_idx)
	{
	  uint64 prefix = my_prefixes[prefix_idx];
	  uint64 target_prefix = prefix ^ s1mask;
	  int target_mpi_rank =
	    get_mpi_rank_spinhalf(target_prefix, mpi_size_);

	    
	  uint64 prefix_size = my_prefix_size[prefix_idx];
	  uint64 prefix_offset = my_prefix_size[prefix_offset];

	  int prefix_n_upspins = popcnt(prefix);
	  int postfix_n_upspins = n_upspins - prefix_n_upspins;
	  assert(postfix_n_upspins >= 0);

	  LinTable& lintable = lintables_[n_sites][postfix_n_upspins];
	  assert(lintable.size() == prefix_size);

	  uint64 begin = prefix_offset;
	  uint64 end = prefix_offset + prefix_size;
	  for (uint64 idx=begin; idx<end; ++idx)
	    {
	      uint64 postfix = lintable.state(idx);
	      if (popcnt(postfix & flipmask)==1)
		{
		  uint64 postfix_flipped = postfix ^ flipmask;
		  uint64 new_idx = begin + lintable.index(postfix_flipped);
		  y[new_idx] += Jmult * x[idx]
		}
	    }
	  
	}  // loop over prefix indices   
    }  // loop over s1, s2, J for postfixes

  
}
