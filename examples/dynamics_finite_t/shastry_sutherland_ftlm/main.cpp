#include <extern/clara/clara.hpp>
#include <xdiag/all.h>

void parse_cmdline(std::string &latfile, double &J, double &Jd, int &nup,
                   std::string &k, std::string &q, int &seed, int &niter,
                   std::string &outfile, std::string &dumpdir, int &argc,
                   char **argv) {

  bool showhelp = false;
  auto parser =
      clara::Opt(latfile,
                 "latfile")["-l"]["--latfile"]("name of lattice file") |
      clara::Opt(J, "J")["-j"]["--J"]("J") |
      clara::Opt(Jd, "Jd")["-d"]["--Jd"]("Jd") |
      clara::Opt(nup, "nup")["-n"]["--nup"]("number of upspins") |
      clara::Opt(k, "k")["-k"]["--k"]("space group irrep") |
      clara::Opt(q, "q")["-q"]["--q"]("momentum q of spin operator") |
      clara::Opt(seed, "seed")["-s"]["--seed"]("random seed") |
      clara::Opt(niter,
                 "niter")["-i"]["--niter"]("number of Lanczos iterations") |
      clara::Opt(outfile, "outfile")["-o"]["--outfile"]("name of outfile") |
      clara::Opt(dumpdir, "dumpdir")["-v"]["--dumpdir"]("name of dumpdir") |
      clara::Help(showhelp);

  auto cmd_args = parser.parse(clara::Args(argc, argv));
  if (!cmd_args) {
    std::cerr << "Error in command line: " << cmd_args.errorMessage()
              << std::endl;
    exit(EXIT_FAILURE);
  } else if (showhelp) {
    parser.writeToStream(std::cout);
    exit(EXIT_SUCCESS);
  } else {
    std::cout << "latfile: " << latfile << std::endl
              << "J      : " << J << std::endl
              << "Jd     : " << Jd << std::endl
              << "nup    : " << nup << std::endl
              << "k      : " << k << std::endl
              << "q      : " << q << std::endl
              << "seed   : " << seed << std::endl
              << "niter  : " << niter << std::endl
              << "outfile: " << outfile << std::endl
              << "dumpdir: " << dumpdir << std::endl
              << "-----" << std::endl;
  }
}

int main(int argc, char **argv) {
  using namespace xdiag;
  using namespace arma;
  using fmt::format;

  say_hello();

  // Parse input arguments
  std::string latfile;
  double J = 1;
  double Jd = 1;
  int nup = 0;
  std::string k = "0";
  std::string q = "0";
  int seed = 1;
  int niter = 200;
  std::string outfile;
  std::string dumpdir;
  parse_cmdline(latfile, J, Jd, nup, k, q, seed, niter, outfile, dumpdir, argc,
                argv);

  auto lfile = FileToml(latfile, 'r');
  auto ofile = FileH5(outfile, "w!");

  auto bonds = BondList(lfile["Interactions"]);
  bonds["J"] = J;
  bonds["Jd"] = Jd;
  auto group = PermutationGroup(lfile["Symmetries"]);
  auto irrep_k = Representation(lfile[k]);
  int n_sites = bonds.n_sites();

  Log("Creating block nup: {}, k: {} ...", nup, k);
  auto block_k = Spinhalf(n_sites, nup, group, irrep_k);
  ofile["DimK"] = block_k.size();
  if (block_k.size() == 0){
    ofile["AlphasT"] = vec(zeros(0));
    ofile["BetasT"] = vec(zeros(0));
    ofile["AlphasS"] = vec(zeros(0));
    ofile["BetasS"] = vec(zeros(0));
    ofile["A"] = mat(zeros(0, 0));
    exit(EXIT_SUCCESS);
  }
  
  Log("Creating Lanczos vector matrix V from |r> ...");
  auto rstate = random_state_cplx(block_k, seed);
  auto dump_V = [dumpdir](int iteration, cx_vec const &vec) {
    std::string filename = format("{}/V{}.arm", dumpdir, iteration);
    tic();
    Log("writing V vector {} -> {}", iteration, filename);
    vec.save(filename);
    toc();
  };
  auto T = lanczos_vector_apply_inplace(bonds, rstate, dump_V, niter);
  ofile["AlphasT"] = T.alphas();
  ofile["BetasT"] = T.betas();

  Log("Creating block nup: {}, k+q: {}+{}", nup, k, q);
  auto irrep_q = Representation(lfile[q]);
  auto irrep_k_q = irrep_k * irrep_q;
  auto block_k_q = Spinhalf(n_sites, nup, group, irrep_k_q);
  ofile["DimKQ"] = block_k_q.size();

  Log("Computing ground state |gs> and S(q)|gs>...");
  auto gs = groundstate(bonds, block_k);
  mat coords = lfile["Coordinates"].as<mat>();
  auto qq = lfile[q + std::string(".momentum")].as<vec>();
  BondList S_of_q;
  for (int site = 0; site < n_sites; ++site) {
    complex phase =
        std::exp(std::complex<double>(0, 1.0) *
                 (qq(0) * coords(site, 0) + qq(1) * coords(site, 1)));
    S_of_q << Bond("SZ", phase / n_sites, site);
  }
  auto s_of_q_gs = StateCplx(block_k_q);
  apply(S_of_q, gs, s_of_q_gs);
  double nrm = norm(s_of_q_gs);
  ofile["SofQNorm"] = nrm;

  if (nrm < 1e-12) {
    Log("Zero norm of S(q)|gs>");
    ofile["AlphasS"] = vec(zeros(0));
    ofile["BetasS"] = vec(zeros(0));
    ofile["A"] = mat(zeros(0, T.size()));
  } else {
    Log("Creating Lanczos vector matrix W from S(q)|g.s.>...");
    auto dump_W = [dumpdir, q](int iteration, cx_vec const &vec) {
      std::string filename = format("{}/W{}.arm", dumpdir, iteration);
      Log("writing W vector {} (q={}) -> {}", iteration, q, filename);
      vec.save(filename);
    };
    auto S = lanczos_vector_apply_inplace(bonds, s_of_q_gs, dump_W, niter);
    ofile["AlphasS"] = S.alphas();
    ofile["BetasS"] = S.betas();

    auto A = cx_mat(S.size(), T.size());
    for (int n = 0; n < T.size(); ++n) {
      Log("Computing W S(q) Vdag (q={}) -> n={}", q, n);
      tic();
      auto v = State(block_k);
      v.vector().load(format("{}/V{}.arm", dumpdir, n));
      auto S_of_q_v = State(block_k_q);
      apply(S_of_q, v, S_of_q_v);

      for (int m = 0; m < S.size(); ++m) {
        auto w = State(block_k_q);
        w.vector().load(format("{}/W{}.arm", dumpdir, m));
        A(m, n) = dot(w, S_of_q_v);
      }
      toc();
    }
    ofile["A"] = A;
  }

  return EXIT_SUCCESS;
}
