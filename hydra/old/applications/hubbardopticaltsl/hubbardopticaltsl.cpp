#include <cstdlib>
#include <vector>
#include <utility>
#include <fstream>
#include <iomanip>
#include <unistd.h>

#include <mpi.h>

#include <lila/allmpi.h>
#include <hydra/allmpi.h>

#include "hubbardopticaltsl.options.h"

#include <iomanip>
#include <locale>

template<class T>
std::string FormatWithCommas(T value)
{
    std::stringstream ss;
    ss.imbue(std::locale(""));
    ss << std::fixed << value;
    return ss.str();
}

int main(int argc, char* argv[])
{
  using hydra::hilbertspaces::hubbard_qn;
  using hydra::hilbertspaces::Hubbard;
  using hydra::models::HubbardModelMPI;
  using namespace hydra::operators;
  using hydra::dynamics::continued_fraction;
  using namespace lila;

  MPI_Init(&argc, &argv); 
  int mpi_rank, mpi_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  // Parse input
  std::string outfile;
  std::string latticefile;
  std::string couplingfile;
  std::string corrfile;
  int nup = -1;
  int ndown = -1;
  int iters = 100;
  int verbosity = 1;
  int seed = 1;

  parse_cmdline(outfile, latticefile, couplingfile, corrfile, nup, ndown, iters, verbosity, seed, argc, argv);


  if ((verbosity >= 1) && (mpi_rank == 0))  
    printf("Using %d MPI tasks\n", mpi_size);

  // Open outfile
  std::ofstream of;
  if (outfile != "")
    {
      of.open(outfile, std::ios::out);
      if(of.fail()) 
	{
	  if (mpi_rank == 0)
	    {
	      std::cerr << "Error in opening outfile: " 
			<< "Could not open file with filename ["
			<< outfile << "] given. Abort." << std::endl;
	      MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	}
    }

  /////////////////////////////////////////////////
  // Create Hamiltonian
  BondList bondlist = read_bondlist(latticefile);
  BondList hopping_list = bondlist.bonds_of_type("HUBBARDHOP");
  BondList interaction_list = bondlist.bonds_of_type("HUBBARDV");
  if (couplingfile == "") couplingfile = latticefile;
  Couplings couplings = read_couplings(couplingfile);

  if ((verbosity >= 1) && (mpi_rank == 0))
    {
      for (auto bond : hopping_list)
	printf("hopping %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto bond : interaction_list)
	printf("interaction %s %d %d\n", bond.coupling().c_str(), 
	       bond.sites()[0], bond.sites()[1]);
      for (auto c : couplings)
	printf("coupling %s %f %fj\n", c.first.c_str(), std::real(c.second), std::imag(c.second));
    }

  int n_sites = bondlist.n_sites();
  hubbard_qn qn;
  if ((nup == -1)  || (ndown == -1))
    qn = {n_sites/2, n_sites/2};
  else qn = {nup, ndown};

  // Create infrastructure for Hubbard model
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Creating Hubbard model for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  double t1 = MPI_Wtime();
  auto model = HubbardModelMPI<complex,uint32>(bondlist, couplings, qn);
  double t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time init: %3.4f\n", t2-t1); 

  // Define multiplication function
  auto H = [&model, &mpi_rank, &verbosity](const VectorMPI<complex>& v, 
							      VectorMPI<complex>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    model.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    // if ((verbosity >= 2) && (mpi_rank == 0))
    //   printf("iter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };
  //
  /////////////////////////////////////////////////


  /////////////////////////////////////////////////
  // Create current operator
  std::vector<std::pair<int, int>> correlation_list;
  if (corrfile == "") 
    {
      if (mpi_rank == 0)
	    {
	      std::cerr << "Error in opening corr: " 
			<< "Could not open file with filename ["
			<< corrfile << "] given. Abort." << std::endl;
	      MPI_Abort(MPI_COMM_WORLD, 3);
	    }
    }  
  BondList curr_bondlist = read_bondlist(corrfile);
  
  // Create infrastructure for current operator
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("Creating current operator for n_upspins=%d, n_downspins=%d...\n",
	   qn.n_upspins, qn.n_downspins);
  t1 = MPI_Wtime();
  for (auto bond : curr_bondlist)
    if ((verbosity >= 1) && (mpi_rank == 0))
      printf("curr %s %d %d\n", bond.coupling().c_str(), 
	     bond.sites()[0], bond.sites()[1]);

  Couplings curr_couplings;
  curr_couplings["C"] = 1;
  auto current = HubbardModelMPI<complex,uint32>(curr_bondlist, curr_couplings, qn);
  t2 = MPI_Wtime();
  if ((verbosity >= 1) && (mpi_rank == 0))
    printf("time curr init: %3.4f\n", t2-t1); 
  // Define multiplication function
  auto A = [&current, &mpi_rank, &verbosity](const VectorMPI<complex>& v, 
					     VectorMPI<complex>& w) {
    static int iter=0;
    double t1 = MPI_Wtime();
    bool verbose = ((iter == 0) && verbosity >=1);
    current.apply_hamiltonian(v, w, verbose);
    double t2 = MPI_Wtime();
    // if ((verbosity >= 2) && (mpi_rank == 0))
    //   printf("iter: %d, time MVM: %3.4f\n", iter, t2-t1); 
    ++iter;
  };
  //
  //////////////////////////////////////////////////


  //////////////////////////////////////////////////
  // Create random start state
  int random_seed = seed + 1234567*mpi_rank;
  uint64 local_dim = model.local_dim();
  VectorMPI<complex> startstate(local_dim);
  normal_dist_t<complex> dist(0., 1.);
  normal_gen_t<complex> gen(dist, random_seed);
  Random(startstate, gen, true);
  Normalize(startstate);
  //
  //////////////////////////////////////////////////

  auto zeros = startstate;
  Zeros(zeros);

  std::vector<VectorMPI<complex>> V;
  std::vector<VectorMPI<complex>> W;

  V.push_back(startstate);
  W.push_back(startstate);
  
  // Two-sided Lanczos
  for (int i=0; i<iters; ++i)
    {
      
      /////////////////////////////////////////////
      auto v = zeros;
      A(V[V.size()-1], v);

      // orthonormalize
      for (int j=0; j<(int)V.size(); ++j)
      	{
      	  auto bv = Dot(V[j], v);
	  v -= bv * V[j];
	}
      auto nv = Norm(v);

      v /= (complex)nv;
      V.push_back(v);

      v = zeros;
      H(V[V.size()-2], v);

      // orthonormalize
      for (int j=0; j<(int)V.size(); ++j)
      	{
      	  auto bv = Dot(V[j], v);
	  v -= bv * V[j];
	}
      nv = Norm(v);
      v /= (complex)nv;
      V.push_back(v);
      /////////////////////////////////////////////
    


      /////////////////////////////////////////////
      auto w = zeros;
      H(W[W.size()-1], w);

      // orthonormalize
      for (int j=0; j<(int)W.size(); ++j)
      	{
      	  auto bw = Dot(W[j], w);
	  w -= bw * W[j];
	}
      auto nw = Norm(w);

      w /= (complex)nw;
      W.push_back(w);

      w = zeros;
      A(W[W.size()-2], w);

      // orthonormalize
      for (int j=0; j<(int)W.size(); ++j)
      	{
      	  auto bw = Dot(W[j], w);
	  w -= bw * W[j];
	}
      nw = Norm(w);
      w /= (complex)nw;
      W.push_back(w);
      /////////////////////////////////////////////




      // // orthogonalize against previous vecs
      // auto WGS = W;
      // WGS = lila::gramschmidt(WGS);
      // auto VGS = V;
      // VGS = lila::gramschmidt(VGS);

      // for (int j=0; j<=i; ++j)
      // 	{

      // 	  // Two-sided Lanczos
      // 	  auto nw = Norm(WGS[j]);
      // 	  auto nv = Norm(VGS[j]);
      // 	  auto bv = Dot(WGS[j], v) / (nw * nw);
      // 	  auto bw = Dot(VGS[j], w) / (nv * nv);
      // 	  v -= bv * WGS[j];
      // 	  w -= bw * VGS[j];

      // 	  // if (mpi_rank==0) printf("i: %d, j %d, bv: (%f, %f), bw: (%f, %f)\n", 
      // 	  // 			  i, j, std::real(bv), std::imag(bv), std::real(bw), std::imag(bw));
      // 	  // if (mpi_rank==0) printf("i: %d, j %d, bv: (%.4f, %.4f), bw: (%.4f, %.4f)\n", 
      // 	  // 			  i, j, std::real(bv), std::imag(bv), std::real(bw), std::imag(bw));
      // 	}

      // auto WGS = W;
      // WGS.insert(WGS.end(), V.begin(), V.end());
      // for (int j=0; j<(int)WGS.size(); ++j)
      // 	{

      // 	  auto nw = Norm(WGS[j]);
      // 	  auto nv = Norm(WGS[j]);
      // 	  auto bv = Dot(WGS[j], v) / (nw * nw);
      // 	  auto bw = Dot(WGS[j], w) / (nv * nv);
      // 	  v -= bv * WGS[j];
      // 	  w -= bw * WGS[j];
      // 	}


  
      // // normalize
      // auto delta = Dot(v, w);
      // auto phase = delta / std::abs(delta);
      // auto nv = sqrt(std::abs(delta));
      // auto nw = sqrt(std::abs(delta)) * phase;

      // auto nv = Norm(v);
      // auto nw = Norm(w);

      // if (mpi_rank==0) printf("i: %d, nv: %.4f, nw: %.4f\n", i, nv, nw);
      
      // v /= (complex)nv;
      // w /= (complex)nw;

      // // for (int j=0; j<=i; ++j)
      // // 	{
      // // 	  auto d = Dot(v, W[j]);
      // // 	  if (mpi_rank==0)
      // // 	    LilaPrint(d);
      // // 	}
      // // if (mpi_rank==0) printf("\n\n\n");
  

      // V.push_back(v);
      // W.push_back(w);			 
    }

  // Compute matrices
  int size = (int)V.size();
  auto T = Matrix<complex>(size, size);
  auto S = Matrix<complex>(size, size);
  auto O = Matrix<complex>(size, size);

  auto TV = Matrix<complex>(size, size);
  auto SW = Matrix<complex>(size, size);
  auto OV = Matrix<complex>(size, size);
  auto OW = Matrix<complex>(size, size);

  // for (int i=0; i<iters; ++i)
  //   for (int j=0; j<iters; ++j)
  //     {
  // 	auto tmp = zeros;
  // 	H(V[j], tmp);
  // 	T(i,j) = Dot(W[i], tmp);
  // 	tmp = zeros;
  // 	A(W[j], tmp);
  // 	S(i,j) = Dot(V[i], tmp);
  // 	O(i,j) = Dot(V[i], W[j]);


  // 	OV(i,j) = Dot(V[i], V[j]);
  // 	OW(i,j) = Dot(W[i], W[j]);
  // 	tmp = zeros;
  // 	H(V[j], tmp);
  // 	TV(i,j) = Dot(V[i], tmp);
  // 	tmp = zeros;
  // 	A(W[j], tmp);
  // 	SW(i,j) = Dot(W[i], tmp);
  //     }
  // if (mpi_rank==0)
  //   {
  //     LilaPrint(O);
  //     LilaPrint(OV);
  //     LilaPrint(OW);
  //     // LilaPrint(T);
  //     LilaPrint(S);

  //     LilaPrint(TV);
  //     LilaPrint(SW);
  //   }
  printf("here\n");
  for (int i=0; i<(int)V.size(); ++i)
    for (int j=0; j<(int)V.size(); ++j)
      {
	auto tmp = zeros;
	H(W[j], tmp);
	auto d = Dot(V[i], tmp); 
	T(i,j) = std::abs(d)<1e-12 ? 0. : d;

	tmp = zeros;
	A(V[j], tmp);
	d = Dot(W[i], tmp); 
	S(i,j) = std::abs(d)<1e-12 ? 0. : d;
	
	d = Dot(V[i], W[j]);
	O(i,j) = std::abs(d)<1e-12 ? 0. : d;
      }


  if (mpi_rank==0)
    {
      LilaPrint(O);
      LilaPrint(T);
      LilaPrint(S);
    }


  MPI_Finalize();
  return EXIT_SUCCESS;
}
