#pragma once

#include <string>
#include <vector>

namespace xdiag::utils {

template <class T>
void check_if_contained_in(T const &elem, std::vector<T> const &vec,
                           std::string name = "element");

void check_if_file_exists(std::string filename);
void check_if_files_exists(std::vector<std::string> filenames);

} // namespace xdiag::utils
