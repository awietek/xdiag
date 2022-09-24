using ITensors

let
  N = 14

  sites = siteinds("S=1/2",N)
  os = OpSum()
  for j=1:N-1
    os += 1.0,"Sz",j,"Sz",j+1
  end
  for j=1:N
    os += 1.0,"Sx",j;
  end
  H = MPO(os,sites)
  sweeps = Sweeps(30)
  maxdim!(sweeps,10,10,10,20,20,40,80,100,200,200,1000)
  cutoff!(sweeps,1E-16)
  # noise!(sweeps,1E-6)
  psi0_init = randomMPS(sites,linkdims=2)
  energy0,psi0 = dmrg(H,psi0_init,sweeps)
  print(energy0)
end
