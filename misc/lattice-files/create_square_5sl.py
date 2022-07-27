#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import quantipy.lattice
import quantipy.models.genericmodel as gm
from quantipy.utils.geometryutils import _map_to_simtorus

##################################################
# Script to create model files for Square lattices
# for the EDFSL engine. Number of sites can only 
# be a multiple of 4
# prefix of filename is "sqaure.$n_sites"
# postfix is ".lat"
##################################################


# 10 site lattice
lattice = quantipy.lattice.Square(simulation_torus = np.array([[5,0],[0,2]]))
filename = "heisenberg"

#############################
# Create sublattice structure
# DO NOT CHANGE ANYTHING HERE
#############################
sl1 = []
sl2 = []
sl3 = []
sl4 = []
sl5 = []
span = np.array([[5.0, 0.0], [0.0, 1.0]])
for idx, coord in enumerate(lattice.coordinates):

    # sublattice structure with C4 invariance
    sl_vec = _map_to_simtorus([coord], span)[0]
    if np.allclose(sl_vec, np.array([0,0])):
        sl1.append(idx)
    elif np.allclose(sl_vec, np.array([1,0])):
        sl2.append(idx)
    elif np.allclose(sl_vec, np.array([2,0])):
        sl3.append(idx)
    elif np.allclose(sl_vec, np.array([3,0])):
        sl4.append(idx)
    elif np.allclose(sl_vec, np.array([4,0])):
        sl5.append(idx)
  
  
sublattice_structure = np.array([sl1,sl2,sl3,sl4,sl5])
print(sublattice_structure)

#Define model
maxpg = ["C1", "C2"]
model = gm.GenericModel(lattice, maxpg=maxpg)


model.set_nb_interaction(1, int_type='HB', coupling='J')
# model.plot()
# plt.show()


# Create modelfile
model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".5sl.lat", 
            sublattice_structure= sublattice_structure)
