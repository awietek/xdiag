// Copyright 2018 Alexander Wietek - All Rights Reserved.
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

#ifndef HYDRA_OPERATORS_COUPLINGS_
#define HYDRA_OPERATORS_COUPLINGS_

#include <map>
#include <string>
#include <vector>

#include <lila/all.h>
#include <hydra/utils/typedefs.h>

namespace hydra { namespace operators {

    class Couplings
    {
      using iterator_t = typename std::map<std::string, complex>::iterator;
      using const_iterator_t = typename std::map<std::string, complex>::const_iterator;
    public:

      Couplings() = default;

      /// Create Couplings from map from string to complex
      explicit Couplings(const std::map<std::string, complex>& couplings)
	: couplings_(couplings)
      { }

      /// returns vector with names of couplings
      std::vector<std::string> couplings() const;

      /// Accessor to coupling value
      complex operator[] (const std::string& name) const 
      { return couplings_.find(name)->second; }
      complex& operator[] (const std::string& name) 
      { return couplings_[name]; }

      /// check if a given coupling name is defined
      bool defined (const std::string& name) const 
      { return couplings_.find(name) != couplings_.end(); }

      /// checks whether coupling has zero imaginary part    
      bool is_real(const std::string& name) const 
      { return lila::close(std::abs(std::imag(couplings_.find(name)->second)), 0.); }

      /// checks whether all couplings have zero imaginary part    
      bool all_real() const 
      {
	for(auto name_val : couplings_)
	  {
	    if (!is_real(name_val.first))
	      return false;
	  }
	return true;	
      }

      /// returns real part of coupling 
      double real(const std::string& name) const 
      { return std::real(couplings_.find(name)->second); }

      /// returns imaginary part of coupling 
      double imag(const std::string& name) const 
      { return std::imag(couplings_.find(name)->second); }


      iterator_t begin() { return couplings_.begin(); }
      iterator_t end() { return couplings_.end(); }
      const_iterator_t begin() const { return couplings_.begin(); }
      const_iterator_t end() const { return couplings_.end(); }
      const_iterator_t cbegin() const { return couplings_.cbegin(); }
      const_iterator_t cend() const { return couplings_.cend(); }

    private:
      std::map<std::string, complex> couplings_;
    };


    /*!
      Reads Couplings from file from filename

      Format should be 
      [Couplings]
      J1 1.1
      J2 2.2
      ...
      reads a block starting with [Couplings] or [couplings] until next 
      [...] block or if [...] is not found until end of file
    */
    Couplings read_couplings(std::string filename);

   
      
  }  // namespace operators
}  // namespace hydra

#endif
