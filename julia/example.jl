using YxDiag

# block = Spinhalf(2, 1)
# bond = Bond("HB", "J", [1, 2])
# op = BondList([bond], Dict("J"=> 1.0))

# e0, gs = eig0(op, block)

# @show e0
# @show col(gs)
# make_complex!(gs)
# @show col(gs)
# @show col(real(gs))
# @show col(imag(gs))


N=4
t=1
U=0

block = Electron(N, N÷2, N÷2)
bonds = []
for i in 1:N
    push!(bonds, Bond("HOP", "t", [i, mod1(i+1, N)]))
end
@show bonds
H = BondList(bonds, Dict("t"=> t, "U"=>U))
e0, gs = eig0(H, block)

psi0 = ProductState(["Up", "Dn", "Up", "Dn"])
@show e0
