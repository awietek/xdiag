using Pkg; Pkg.instantiate()
using XDiag
using HDF5

function main()
  say_hello();
  filename = "../../../misc/data/examples_output/spinhalf_chain_structure_factor.h5"
  outfile = h5open(filename, "w")

  # define model
  N = 16
  nup = div(N,2)
  ops = OpSum()
  for i in 1:N
    ops += "J" * Op("SdotS", [i, mod1(i + 1, N)]);
  end
  ops["J"] = 1.0
  @show ops

  # Create the permutation group
  @show perm = Permutation([mod1(s+1,N) for s in 1:N])
  @show group = PermutationGroup([perm ^ s for s in 0:N-1])

  # Create the irreps at momenta k
  irreps = Representation[]
  for k in 0:N-1
    phase = exp(2im * pi * k / N)
    characters = [phase ^ s for s in 0:N-1]
    irrep = Representation(group, characters)
    # @show irrep
    push!(irreps,irrep)
  end

  # compute groundstate (known to be at k=0)
  println("Computing ground state ...")
  block = Spinhalf(N, nup, irreps[1])
  e0, gs = eig0(ops, block)
  make_complex!(gs)
  println("done.")
  println("Ground state energy: $e0")
  outfile["e0"] = e0;

  # loop through momenta
  for k in 1:N
    println("coumputing structure factor at q=$k")
    S_q = symmetrize(Op("Sz", 1), irreps[k])
    Av = apply(S_q, gs)
    nrm = norm(Av)
    Av /= nrm

    res = eigvals_lanczos_inplace(ops, Av)
    outfile["$k/norm"] = nrm
    outfile["$k/alphas"] = res.alphas
    outfile["$k/betas"] = res.betas
    outfile["$k/eigs"] = res.eigenvalues
  end
end

main()