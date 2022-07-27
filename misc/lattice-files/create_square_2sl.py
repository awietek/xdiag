#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import quantipy.lattice
import quantipy.models.genericmodel as gm
from quantipy.utils.geometryutils import _map_to_simtorus

lattice = quantipy.lattice.Square(simulation_torus = np.array([[2,2],[2,-2]]))
filename = "heisenberg"

#############################
# Create sublattice structure
# DO NOT CHANGE ANYTHING HERE
#############################
sl1 = []
sl2 = []
span = np.array([[2.0, 0.0], [1.0, 1.0]])
# span = np.array([[2.0, 0.0], [0, 1.0]])
for idx, coord in enumerate(lattice.coordinates):
    print idx, coord , int(coord[0] + 2 *coord[1]) % 4, _map_to_simtorus([coord], span)[0]
    
    # sublattice structure with C4 invariance
    sl_vec = _map_to_simtorus([coord], span)[0]
    if np.allclose(sl_vec, np.array([0,0])):
        sl1.append(idx)
    else:
        sl2.append(idx)
 

  
  
sublattice_structure = np.array([sl1,sl2])
print(sublattice_structure)

model = gm.GenericModel(lattice)


model.set_nb_interaction(1, int_type='HB', coupling='J')

model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".2sl.lat", 
            sublattice_structure= sublattice_structure)

