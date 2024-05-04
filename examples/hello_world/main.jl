using XDiag
say_hello()
block = Electron(2, 1, 1);
bond = Bond("HB", 1.0, [0]);
mat = matrix(bond, block);
