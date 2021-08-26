#include "couplings.h"
#include <fstream>
#include <iostream>

namespace hydra {

std::vector<std::string> Couplings::couplings() const {
  std::vector<std::string> cps;
  for (auto c : couplings_)
    cps.push_back(c.first);
  return cps;
}

Couplings read_couplings(std::string filename) {

  std::map<std::string, complex> coupling_map;

  // Open file and handle error
  std::ifstream File(filename.c_str());
  if (File.fail()) {
    std::cerr << "Error in read_couplings: "
              << "Could not open file with filename [" << filename
              << "] given. Abort." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Advance to coupling lines
  std::string tobeparsed;
  getline(File, tobeparsed);
  while ((tobeparsed.find("[Couplings]") == std::string::npos) &&
         (tobeparsed.find("[couplings]") == std::string::npos))
    getline(File, tobeparsed);

  // read lines until '[' is found or else until EOF
  while (std::getline(File, tobeparsed)) {
    if ((tobeparsed.find('[') != std::string::npos))
      break;
    // Parse line
    std::string name;
    complex val;
    std::stringstream stream(tobeparsed);
    stream >> name;
    stream >> val;
    coupling_map[name] = val;
    if (!File.good())
      break;
  }

  return Couplings(coupling_map);
}

bool is_complex(Couplings const &cpls) { return !cpls.all_real(); }
bool is_real(Couplings const &cpls) { return !is_complex(cpls); }

} // namespace hydra
