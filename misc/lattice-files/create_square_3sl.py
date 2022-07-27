#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import quantipy.lattice
import quantipy.models.genericmodel as gm
from quantipy.utils.geometryutils import _map_to_simtorus



# 9 site lattice 
lattice = quantipy.lattice.Square(simulation_torus = np.array([[3,0],[0,3]]))
filename = "heisenberg" 



sl1 = []
sl2 = []
sl3 = []
span = np.array([[3.0, 0.0], [1.0, 1.0]])
# span = np.array([[2.0, 0.0], [0, 1.0]])
for idx, coord in enumerate(lattice.coordinates):

    # sublattice structure with C4 invariance
    sl_vec = _map_to_simtorus([coord], span)[0]
    if np.allclose(sl_vec, np.array([0,0])):
        sl1.append(idx)
    elif np.allclose(sl_vec, np.array([1,0])):
        sl2.append(idx)
    else:
        sl3.append(idx)
 

  
  
sublattice_structure = np.array([sl1,sl2,sl3])
print(sublattice_structure)

#Define model
maxpg = ["C1", "C2", "C4s", "C4^3s"]
model = gm.GenericModel(lattice, maxpg=maxpg)


model.set_nb_interaction(1, int_type='HB', coupling='J')
# model.plot()
# plt.show()


# Create modelfile
model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".3sl.lat", 
            sublattice_structure= sublattice_structure)

