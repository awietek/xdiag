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

#ifndef HYDRA_SYMMETRIES_SPACEGROUP_
#define HYDRA_SYMMETRIES_SPACEGROUP_

#include <string>
#include <vector>

#include <hydra/symmetries/symmetrydetail.h>
#include <hydra/hilbertspaces/siteoperations.h>

namespace hydra { namespace symmetries {

    /*!
      SpaceGroup is a container for lattice symmetries, with the ability
      to apply a symmetry onto a given state

      Usage:
      @code
      #include <hydra/all.h>
      
      using namespace hydra::symmetries;

      std::vector<std::vector<int>> symmetries;
      symmetries.push_back({0, 1, 2, 3});
      symmetries.push_back({1, 2, 3, 0});
      symmetries.push_back({2, 3, 0, 1});
      symmetries.push_back({3, 0, 1, 2});
  
      SpaceGroup sg(symmetries);
      Print(sg);

      SpaceGroup sub_sg = sg.subgroup({0, 2});
      Print(sub_sg);

      SpaceGroup square_sg = read_spacegroup("misc/lattice-files/square.16.spinlessfermions.pbc");
      Print(square_sg);
      @endcode
    */
    class SpaceGroup
    {
    public:
      /*!
	Constructor of SpaceGroup from a vector of int vectors, containing 
	the lattice permutation

	@param symmetries vector of int vector containing permutations
      
	Usage example:
	@code
	#include <hydra/all.h>
      
	using namespace hydra::symmetries;

	std::vector<std::vector<int>> symmetries;
	symmetries.push_back({0, 1, 2, 3});
	symmetries.push_back({1, 2, 3, 0});
	symmetries.push_back({2, 3, 0, 1});
	symmetries.push_back({3, 0, 1, 2});
	SpaceGroup sg(symmetries);
	@endcode
      */
      SpaceGroup(const std::vector<std::vector<int>>& symmetries);

      template <class hilbertspace_t>
      inline typename hilbertspace_t::state_t 
      apply(const int& n_sym, const typename hilbertspace_t::state_t& state) 
	const
      {
	using hydra::hilbertspaces::get_site_val;
	using hydra::hilbertspaces::set_site_val;
	return detail::apply_permutation
	  (state, n_sites_, symmetries_internal_.data() + n_sym*n_sites_, 
	   get_site_val<hilbertspace_t>, set_site_val<hilbertspace_t>);
      }


      template <class hilbertspace_t>
      inline double fermi_sign
      (const int& n_sym, const typename hilbertspace_t::state_t& state) 
	const
      { 
	return detail::fermi_sign<hilbertspace_t>
	  (state, n_sites_, symmetries_internal_.data() + n_sym*n_sites_);  
      }

      /*!
	Return a subgroup of the space group 

	@param symmetry_numbers numbers of symmetries contained in subgroup
      
	Usage example:
	@code
	using namespace hydra::symmetries;

	std::vector<std::vector<int>> symmetries;
	symmetries.push_back({0, 1, 2, 3});
	symmetries.push_back({1, 2, 3, 0});
	symmetries.push_back({2, 3, 0, 1});
	symmetries.push_back({3, 0, 1, 2});
  
	SpaceGroup sg(symmetries);
	Print(sg);

	SpaceGroup sub_sg = sg.subgroup({0, 2});
	Print(sub_sg);
	@endcode
      */
      SpaceGroup subgroup(const std::vector<int>& symmetry_numbers) const;
      
      int n_sites() const { return n_sites_; }
      int n_symmetries() const { return n_symmetries_; }
      const std::vector<std::vector<int>>& symmetries() const { return symmetries_; }
      
    private:
      int n_sites_;
      int n_symmetries_;
      std::vector<std::vector<int>> symmetries_;
      std::vector<int> symmetries_internal_;  // size = n_symmetries_*n_sites_
    };

    void Print(const SpaceGroup& group);
    SpaceGroup read_spacegroup(std::string filename);

  }
}

#endif
