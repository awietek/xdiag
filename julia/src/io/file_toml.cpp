#include "file_toml.hpp"

namespace xdiag::julia {
void define_file_toml(jlcxx::Module &mod) {
  mod.add_type<FileToml>("cxx_FileToml")
      .constructor<>()
      .constructor<std::string>()
      .method("defined", [](FileToml const &fl, std::string key) {
        JULIA_XDIAG_CALL_RETURN(fl.defined(key));
      });
}

} // namespace xdiag::julia
