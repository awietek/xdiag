using XDiag
using HDF5

function main()
  say_hello()

  # IO
  latticeInput =
       "../../misc/data/square.8.hubbard.ttprime.toml"
  lfile = FileToml(latticeInput)

  filename =
      "../../misc/data/examples_output/hubbard_greens_f.h5"
  outfile = h5open(filename, "w")

  # Define the Hubbard model
  N = 8
  nup = div(N,2)
  ndn = div(N,2)
  ops = read_opsum(lfile, "Interactions")
  ops["Ty"] = 1.
  ops["Tx"] = 1.
  ops["Tprime+"] = 0.3
  ops["Tprime-"] = 0.3
  ops += "U" * Op("HubbardU")
  ops["U"] = -6.
  @show(ops)
  opsTBC = ops

  irrep = read_representation(lfile, "Gamma.C1.A")

  # compute groundstate (known to be at k=0)
  println("Computing ground state ...")
  block = Electron(N, nup, ndn, irrep)
  e0, gs = eig0(ops, block)
  make_complex!(gs)
  println("done.")
  println("Ground state energy: $e0")
  outfile["e0"] = e0

  # loop through momenta
  irreps = Dict{String,Vector{Float64}}(
      "Gamma.C1.A"=> [0., 0.],
      "M.C1.A"=> [3.1415926535897931, 3.1415926535897931],
      "Sigma0.C1.A"=> [1.5707963267948966, 1.5707963267948966],
      "Sigma1.C1.A"=> [1.5707963267948966, -1.5707963267948966],
      "Sigma2.C1.A"=> [-1.5707963267948966, 1.5707963267948966],
      "Sigma3.C1.A"=> [-1.5707963267948966, -1.5707963267948966],
      "X0.C1.A"=> [3.1415926535897931, 0.0000000000000000],
      "X1.C1.A"=> [0.0000000000000000, 3.1415926535897931])
  for (name, momentum) in irreps 
    println("considering irrep $name, momentum $momentum")
    aq_irrep = read_representation(lfile, name)
    c_q = symmetrize(Op("Cup", 1), aq_irrep)
    Av = apply(c_q, gs)
    @show nrm = norm(Av)
    Av /= nrm

    res = eigvals_lanczos_inplace(ops, Av)
    @show res.eigenvalues[1]
    outfile["$name/norm"] = nrm
    outfile["$name/alphas"] = res.alphas
    outfile["$name/betas"] = res.betas
    outfile["$name/eigs"] = res.eigenvalues

    # compare to twisted boundary condition calculation,
    # cf. Zemljic & Prelovsek, PRB 75, 104514 (2007)
    # Tohyama, PRB 70, 174517 (2004).
    opsTBC["Tx"] = exp(1im * momentum[1])
    opsTBC["Ty"] = exp(1im * momentum[2])
    opsTBC["Tprime+"] =
        0.3 * exp(1im * (momentum[1] + momentum[2])) # r_ij = (1,1)
    opsTBC["Tprime-"] =
        0.3 * exp(1im * (momentum[1] - momentum[2])) # r_ij = (1,-1)
    c_q_tbc = symmetrize(Op("Cup", 1), irrep)
    Av = apply(c_q_tbc, gs)
    @show nrm = norm(Av)
    Av /= nrm

    res_tbc = eigvals_lanczos_inplace(opsTBC, Av)
    @show res_tbc.eigenvalues[1]
  end
end

main()
