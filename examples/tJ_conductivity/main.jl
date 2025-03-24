using XDiag
using HDF5

function main()
  say_hello()

  # IO
  latticeInput = "../../misc/data/square.16.tJ.toml"
  lfile = FileToml(latticeInput)

  filename = "../../misc/data/examples_output/tJ_conductivity_jll.h5"
  outfile = h5open(filename, "w")

  # Lanczos parameters
  precision = 0.0 # turns off precision checking
  maxiters = 200

  # Define the model
  N = 16
  nup = 8
  ndn = 7 # one hole
  ops = read_opsum(lfile, "Interactions")
  t = 1.0
  J = 0.3
  ops["T"] = t
  ops["J"] = J
  @show(ops)
  irrep = read_representation(lfile, "Gamma.C1.A")

  # compute groundstate
  println("Computing ground state ...")
  block = tJ(N, nup, ndn, irrep)
  e0, gs = eig0(ops, block)
  println("done.")
  println("Ground state energy: $e0")
  outfile["e0"] = e0

  current = symmetrize(1.0 * im * Op("Hop", [1, 2]), irrep)
  Av = apply(current, gs)
  @show nrm = norm(Av)
  Av /= nrm

  res = eigvals_lanczos_inplace(ops, Av; neigvals=1, precision=precision, max_iterations=maxiters)
  @show res.eigenvalues[1]
  outfile["norm"] = nrm
  outfile["alphas"] = res.alphas
  outfile["betas"] = res.betas
  outfile["eigs"] = res.eigenvalues
end

main()
