#include <xdiag/all.hpp>

using namespace xdiag;

void push_to_ps(ProductState &ps1, ProductState &psA, ProductState &psB);
void Create_RDM(arma::cx_mat &RDM, arma::vec &vstate, Spinhalf &Block,
                Spinhalf &BlockA, Spinhalf &BlockB);
double getVnE(arma::cx_mat &RDM);

int main() try {
  int N = 16;         // chain length
  double Delta = 1.0; // anisotropy

  auto ops = OpSum();
  auto block = Spinhalf(N); // Spin1/2 Block full

  for (int i = 0; i < N; i++) {
    ops += Op("Exchange", {i, (i + 1) % N});
    ops += Delta * Op("SzSz", {i, (i + 1) % N});
  }

  auto [e0, gs] = eig0(ops, block); // get GS with Lanczos

  arma::vec vstate = vector(gs);

  arma::mat vne = arma::mat(
      N / 2, 2,
      arma::fill::zeros); // Array to store EE as a function of subsystem size

  for (int na = 0; na < N / 2;
       na++) // Compute the Entanglement entropy for different subsystems from
             // l=1 to l=N/2, as the entanglement entropy satisfies: S(l) =
             // S(N=l).
  {
    int nb = N - (na + 1);
    auto blockA = Spinhalf((na + 1)); // Spin1/2 Block full
    auto blockB = Spinhalf(nb);       // Spin1/2 Block full
    arma::cx_mat RDM =
        arma::cx_mat(dim(blockA), dim(blockA), arma::fill::zeros);
    Create_RDM(RDM, vstate, block, blockA, blockB);
    vne(na, 0) = (na + 1);
    vne(na, 1) = getVnE(RDM);
  }

  // Construct the filename
  std::string flstring = "EE_XXX_model.Nsites." + std::to_string(N) + ".Delta" +
                         std::to_string((double)Delta) + ".outfile.h5";

  auto save_fl = FileH5(flstring, "w!");
  save_fl["entanglement_entropy"] = vne;
  return 0;
} catch (Error e) {
  error_trace(e);
}

void push_to_ps(ProductState &ps1, ProductState &psA, ProductState &psB) {
  for (auto sub : psA) {
    ps1.push_back(sub);
  }
  for (auto sub : psB) {
    ps1.push_back(sub);
  }
}

void Create_RDM(arma::cx_mat &RDM, arma::vec &vstate, Spinhalf &Block,
                Spinhalf &BlockA, Spinhalf &BlockB) {
  int i = 0;
  int j = 0;
  for (auto p1 : BlockA) {
    j = 0;
    for (auto p2 : BlockA) {
      for (auto p3 : BlockB) // loop over all configurations in subspace that we
                             // are interested in tracing.
      {
        auto p1aux = ProductState();
        auto p2aux = ProductState();
        push_to_ps(p1aux, p1, p3);
        push_to_ps(p2aux, p2, p3);
        int id1 = index(Block, p1aux);
        int id2 = index(Block, p2aux);
        RDM(j, i) += vstate(id2) * std::conj(vstate(id1));
      }
      j++;
    }
    i++;
  }
}

double getVnE(arma::cx_mat &RDM) {
  // Compute the Entanglement entropy
  double vne = 0.0;
  arma::vec pvals;
  arma::cx_mat eigvec;
  arma::eig_sym(pvals, eigvec, RDM); // FULL ED in the RDM

  for (int p = 0; p < pvals.n_elem; p++) {
    if (pvals(p) > 1e-14)
      vne -= pvals(p) * log(pvals(p));
  }
  return vne;
};
