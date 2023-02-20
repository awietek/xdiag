#include "say_hello.h"

#include <hydra/config.h>
#include <iostream>
#include <sstream>

namespace hydra {

std::string version_string() {
  std::stringstream sstr;
  sstr << "Hydra Version: " << HYDRA_VERSION << "\n";
  sstr << "Git hash     : " << git_hash() << "\n";
  sstr << "Directory    : " << HYDRA_DIRECTORY << "\n";
  sstr << "Hostname     : " << HYDRA_HOSTNAME << "\n";
  sstr << "Compiled by  : " << HYDRA_COMPILEDBY << "\n";
  return sstr.str();
}

void say_hello() {
  std::string hello =
      R"(
.---.  .---.   ____     __  ______     .-------.       ____
|   |  |_ _|   \   \   /  /|    _ `''. |  _ _   \    .'  __ `.
|   |  ( ' )    \  _. /  ' | _ | ) _  \| ( ' )  |   /   '  \  \
|   '-(_{;}_)    _( )_ .'  |( ''_'  ) ||(_ o _) /   |___|  /  |
|      (_,_) ___(_ o _)'   | . (_) `. || (_,_).' __    _.-`   |
| _ _--.   ||   |(_,_)'    |(_    ._) '|  |\ \  |  |.'   _    |
|( ' ) |   ||   `-'  /     |  (_.\.' / |  | \ `'   /|  _( )_  |
(_{;}_)|   | \      /      |       .'  |  |  \    / \ (_ o _) /
'(_,_) '---'  `-..-'       '-----'`    ''-'   `'-'   '.(_,_).'
)";
  std::cout << hello << "\n" << version_string() << "\n";
}

void print_version() { std::cout << version_string() << "\n"; }
} // namespace hydra
