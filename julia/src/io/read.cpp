#include "read.hpp"

namespace xdiag::julia {

void define_read(jlcxx::Module &mod) {

  mod.method("cxx_read_permutation_group", [](FileToml file, std::string tag) {
    JULIA_XDIAG_CALL_RETURN(read_permutation_group(file, tag));
  });
  mod.method("cxx_read_representation", [](FileToml file, std::string irrep_tag,
                                           std::string group_tag) {
    JULIA_XDIAG_CALL_RETURN(read_representation(file, irrep_tag, group_tag));
  });
  mod.method("cxx_read_opsum", [](FileToml file, std::string tag) {
    JULIA_XDIAG_CALL_RETURN(read_opsum(file, tag));
  });
}

} // namespace xdiag::julia
