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

#pragma once

namespace hydra { namespace utils {

  template <typename size_type=int>
  class range
  {
  public:
    constexpr range(size_type e)
      : range(0, e)
    { }

    constexpr range(size_type b, size_type e)
      : current_(b),
        end_(e)
    { }

    size_type const&
    operator*() const { return current_; }

    range& operator++() { 
      ++current_; 
      return *this;
    }

    bool operator!=(const range& other) const
    { return current_ != other.current_; }

    range begin() const { return range(current_, end_); }
    range end() const { return range(end_, end_); }

  private:
    size_type current_;
    size_type end_;
  };


  }  // namespace utils
}  // namespace hydra
