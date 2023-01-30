#include <hydra/all.h>

int main(int argc, char** argv) {
  using namespace hydra;

  assert(argc==2);
  std::string lfile = argv[1];
  BondList bonds = read_bondlist(lfile);

  auto tfile = FileToml(lfile + std::string(".toml"), 'w');
  tfile["Interactions"] = bonds;
  tfile.write();
  return EXIT_SUCCESS;
}
