using XDiag
using HDF5

function main()
  say_hello()

  # IO
  latticeInput = "../../misc/data/square.8.hubbard.ttprime.toml"
  lfile = FileToml(latticeInput)

  filename = "../../misc/data/examples_output/hubbard_greens_f.h5"
  outfile = h5open(filename, "w")

  # Define the Hubbard model
  N = 8
  nup = 3
  ndn = 3
  ops = read_opsum(lfile, "Interactions")
  t = 1.0
  tp = 0.0
  ops["Ty"] = t
  ops["Tx"] = t
  ops["Tprime+"] = tp
  ops["Tprime-"] = tp
  ops += "U" * Op("HubbardU")
  ops["U"] = -6.0
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

  function tbc_mesh(kx, ky)
    # compare to twisted boundary condition calculation,
    # cf. Zemljic & Prelovsek, PRB 75, 104514 (2007)
    # Tohyama, PRB 70, 174517 (2004).
    println("\n considering TBC momentum [$kx,$ky]")
    opsTBC["Tx"] = t * exp(1im * kx)
    opsTBC["Ty"] = t * exp(1im * ky)
    opsTBC["Tprime+"] =
      tp * exp(1im * (kx + ky)) # r_ij = (1,1)
    opsTBC["Tprime-"] =
      tp * exp(1im * (kx - ky)) # r_ij = (1,-1)

    e0, gs_tbc = eig0(opsTBC, block)
    # @show e0
    c_q_tbc = symmetrize(sqrt(N) * Op("Cup", 1), irrep)
    Av = apply(c_q_tbc, gs_tbc)
    @show nrm = norm(Av)
    Av /= nrm

    res_tbc = eigvals_lanczos_inplace(opsTBC, Av)
    @show res_tbc.eigenvalues[1]
    return nrm, res_tbc
  end

  # loop through momenta
  irreps = Dict{String,Vector{Float64}}(
    "Gamma.C1.A" => [0.0, 0.0],
    "M.C1.A" => [3.1415926535897931, 3.1415926535897931],
    "Sigma0.C1.A" => [1.5707963267948966, 1.5707963267948966],
    "Sigma1.C1.A" => [1.5707963267948966, -1.5707963267948966],
    "Sigma2.C1.A" => [-1.5707963267948966, 1.5707963267948966],
    "Sigma3.C1.A" => [-1.5707963267948966, -1.5707963267948966],
    "X0.C1.A" => [3.1415926535897931, 0.0000000000000000],
    "X1.C1.A" => [0.0000000000000000, 3.1415926535897931])
  for (name, momentum) in irreps
    println("\n considering irrep $name, momentum $momentum")
    aq_irrep = read_representation(lfile, name)
    c_q = symmetrize(sqrt(N) * Op("Cup", 1), aq_irrep)
    Av = apply(c_q, gs)
    @show nrm = norm(Av)
    Av /= nrm

    res = eigvals_lanczos_inplace(ops, Av)
    @show res.eigenvalues[1]
    outfile["$name/norm"] = nrm
    outfile["$name/alphas"] = res.alphas
    outfile["$name/betas"] = res.betas
    outfile["$name/eigs"] = res.eigenvalues

    tbc_mesh(momentum[1], momentum[2])
  end

  for (k, kx) in enumerate(0:pi/101:pi)
    nrm, res_tbc = tbc_mesh(kx, kx)

    ik = "$(k-1)" # mimics c++ output
    # group = create_group(outfile, ik)
    # subgroup = create_group(group, ik)
    outfile["$ik/$ik/norm"] = nrm
    outfile["$ik/$ik/alphas"] = res_tbc.alphas
    outfile["$ik/$ik/betas"] = res_tbc.betas
    outfile["$ik/$ik/eigs"] = res_tbc.eigenvalues
  end
end

main()
