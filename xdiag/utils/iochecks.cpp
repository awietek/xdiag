// SPDX-FileCopyrightText: 2025 Alexander Wietek <awietek@pks.mpg.de>
//
// SPDX-License-Identifier: Apache-2.0

#include "iochecks.hpp"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <string>

namespace xdiag {
namespace utils {

template <class T>
void check_if_contained_in(const T &elem, const std::vector<T> &vec,
                           std::string name) {
  assert(vec.size() > 0);
  if (std::find(vec.begin(), vec.end(), elem) == vec.end()) {
    std::cerr << "Unknown " << name << ": " << elem << " (expected ";
    for (int i = 0; i < (int)vec.size() - 1; ++i)
      std::cerr << vec[i] << "/";
    std::cerr << vec[vec.size() - 1] << ")\n";
    exit(EXIT_FAILURE);
  }
}

void check_if_file_exists(std::string filename) {
  std::ifstream f(filename.c_str());
  if (!f.good()) {
    std::cerr << "Could not open file with filename [" << filename
              << "] given. Abort." << std::endl;
    exit(EXIT_FAILURE);
  }
}

void check_if_files_exists(std::vector<std::string> filenames) {
  for (auto filename : filenames)
    check_if_file_exists(filename);
}

template void
check_if_contained_in<std::string>(const std::string &elem,
                                   const std::vector<std::string> &vec,
                                   std::string name);

} // namespace utils
} // namespace xdiag
