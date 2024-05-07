using XDiag
say_hello()

block = Spinhalf(2);
    
bonds = BondList()
bonds += Bond("HB", "J", [0])
bonds["J"] = 1.0;

set_verbosity(2);             # set verbosity for monitoring progress
e0 = eigval0(bonds, block);   # compute ground state energy
