using YxDiag

block = Spinhalf(2, 1)
bond = Bond("HB", "J", [1, 2])
op = BondList([bond], Dict("J"=> 1.0))

e0, gs = eig0(op, block)

@show e0
@show col(gs)
make_complex!(gs)
@show col(gs)
@show col(real(gs))
@show col(imag(gs))

