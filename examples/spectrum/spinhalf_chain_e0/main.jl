using XDiag

let 
    n_sites = 16;
    nup = n_sites ÷ 2;
    block = Spinhalf(n_sites, nup);
    
    # Define the nearest-neighbor Heisenberg model
    bonds = BondList()
    for i in 1:n_sites
        bonds += Bond("HB", "J", [i-1, i % n_sites])
    end
    bonds["J"] = 1.0;

    set_verbosity(2);             # set verbosity for monitoring progress
    e0 = eigval0(bonds, block);   # compute ground state energy

    println("Ground state energy: $e0");
end


