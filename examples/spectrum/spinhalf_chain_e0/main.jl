using XDiag

let 
    N = 16;
    nup = N ÷ 2;
    block = Spinhalf(N, nup);
    
    # Define the nearest-neighbor Heisenberg model
    bonds = BondList()
    for i in 1:N
        bonds += Bond("HB", "J", [i-1, i % N])
    end
    bonds["J"] = 1.0;

    set_verbosity(2);             # set verbosity for monitoring progress
    e0 = eigval0(bonds, block);   # compute ground state energy

    println("Ground state energy: $e0");
end


