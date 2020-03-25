#!/usr/bin/env python
import h5py
import numpy as np

ltype = "chain"
L = 6

filename="alpsraw/spectra.{}.{}.task1.out.h5".format(ltype, L)

spectrum = []    
with h5py.File(filename, "r") as f:

    for key in f["spectrum"]["sectors"].keys():
        sector = f["spectrum"]["sectors"][key]
        spectrum += list(sector["energies"][:])
        print(sector["energies"][:])

spectrum.sort()
np.savetxt("spectrum.{}.{}.txt".format(ltype, L), np.array(spectrum))
