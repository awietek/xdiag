#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import quantipy.lattice
import quantipy.models.genericmodel as gm
from quantipy.utils.geometryutils import _map_to_simtorus


lattice = quantipy.lattice.Square(simulation_torus = np.array([[2,2],[2,-2]]))
filename = "heisenberg"


sl1 = []
sl2 = []
sl3 = []
sl4 = []
span = np.array([[2.0, 0.0], [0.0, 2.0]])
# span = np.array([[2.0, 0.0], [0, 1.0]])
for idx, coord in enumerate(lattice.coordinates):
    print(idx, coord , int(coord[0] + 2 *coord[1]) % 4, _map_to_simtorus([coord], span)[0])
    
    # sublattice structure with C4 invariance
    sl_vec = _map_to_simtorus([coord], span)[0]
    if np.allclose(sl_vec, np.array([0,0])):
        sl1.append(idx)
    elif np.allclose(sl_vec, np.array([1,0])):
        sl2.append(idx)
    elif np.allclose(sl_vec, np.array([0,1])):
        sl3.append(idx)
    elif np.allclose(sl_vec, np.array([1,1])):
        sl4.append(idx)
    else:
        print("hsit")
  
  
sublattice_structure = np.array([sl1,sl2,sl3,sl4])
print(sublattice_structure)

#Define model
maxpg = ["C1", "C2", "C4s", "C4^3s"]
model = gm.GenericModel(lattice)


model.set_nb_interaction(1, int_type='HB', coupling='J')

model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".4sl.lat", 
            sublattice_structure= sublattice_structure)


# irreps = model.symmetries()
# print irreps
# irrepnames = irreps.keys()
# print "(\'" + "\' \'".join(irrepnames) + "\')"

