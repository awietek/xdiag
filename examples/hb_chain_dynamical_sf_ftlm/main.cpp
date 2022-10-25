#include <filesystem>
#include <hydra/all.h>

int main(int argc, char **argv) {
  using namespace hydra;
  using namespace arma;

  assert(argc == 5);
  int n_sites = atoi(argv[1]);
  int n_up = atoi(argv[2]);
  int seed = atoi(argv[3]);
  int n_iters = atoi(argv[4]);

  std::string outdir =
      fmt::format("outfiles/N.{}.nup.{}/seed.{}", n_sites, n_up, seed);
  std::string scratchdir = "/scratch/awietek/tmp";
  std::filesystem::create_directories(outdir);
  std::filesystem::create_directories(scratchdir);

  // Create nearest-neighbor Heisenberg model
  BondList bonds;
  for (int s = 0; s < n_sites; ++s) {
    bonds << Bond("HB", "J", {s, (s + 1) % n_sites});
  }
  bonds["J"] = 1.0;

  Log("Creating Lanczos vector matrix V from |r> ...");
  auto block = Spinhalf(n_sites, n_up);
  auto rstate = random_state_cplx(block, seed);

  auto dump_vector = [scratchdir](int iteration, cx_vec const &vec) {
    std::string filename = fmt::format("{}/v{}.arm", scratchdir, iteration);
    Log("writing {}", filename);
    vec.save(filename);
  };
  auto tmat_V =
      lanczos_vector_apply_inplace(bonds, rstate, dump_vector, n_iters);
  tmat_V.alphas().save(hdf5_name(fmt::format("{}/outfile.h5", outdir), "alphas",
                                 hdf5_opts::append));
  tmat_V.betas().save(hdf5_name(fmt::format("{}/outfile.h5", outdir), "betas",
                                hdf5_opts::append));

  // Reset number of iterations in case deflation took place
  n_iters = tmat_V.size();

  auto gs = groundstate(bonds, block);

  for (int q = 0; q < n_sites; ++q) {
    tic();
    Log("Dynamical Lanczos iterations for q={}", q);

    complex phase = exp(2i * pi * q / n_sites);
    BondList S_of_q;
    for (int s = 0; s < n_sites; ++s) {
      S_of_q << Bond("SZ", pow(phase, s) / n_sites, s);
    }

    // Compute < v_j | S(q) | v_i >
    cx_mat S_of_q_mat(n_iters, n_iters);
    for (int i = 0; i < n_iters; ++i) {
      auto vi = zero_state_cplx(block);
      std::string filename = fmt::format("{}/v{}.arm", scratchdir, i);
      vi.vector().load(filename);
      auto Svi = zero_state_cplx(block);
      apply(S_of_q, vi, Svi);

      for (int j = 0; j < n_iters; ++j) {
        auto vj = zero_state_cplx(block);
        std::string filename = fmt::format("{}/v{}.arm", scratchdir, j);
        vj.vector().load(filename);
        S_of_q_mat(i, j) = dot(vj, Svi);
      }
    }
    S_of_q_mat.save(hdf5_name(fmt::format("{}/outfile.h5", outdir),
                              fmt::format("s_of_q_{}", q), hdf5_opts::append));
    HydraPrint(S_of_q_mat);
    toc();
  }

  std::filesystem::remove_all(fmt::format("{}/Vmatrix", outdir));

  return EXIT_SUCCESS;
}
