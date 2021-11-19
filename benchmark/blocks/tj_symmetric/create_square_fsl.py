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
obc=False
if obc:
    filename = "tJ.fsl.obc" 
else:
    filename = "tJ.fsl.pbc" 


# 9 site lattice
lattice = quantipy.lattice.Square(simulation_torus = np.array([[3,0],[0,3]]))
    
    
# # 16 site lattice
# lattice = quantipy.lattice.Square(simulation_torus = np.array([[4,0],[0,4]]))
    
# # 20 site lattice
# lattice = quantipy.lattice.Square(simulation_torus = np.array([[4,2],[-2,4]]))

# # 25 site lattice
# lattice = quantipy.lattice.Square(simulation_torus = np.array([[5,0],[0,5]]))

# # 32 site lattice
# lattice = quantipy.lattice.Square(simulation_torus = np.array([[4,4],[4,-4]]))

# # 36 site lattice
# lattice = quantipy.lattice.Square(simulation_torus = np.array([[6,0],[0,6]]))


########################
# Print and Plot lattice
########################

for idx, coord in enumerate(lattice.coordinates):
    print(idx, coord)
# perms, irreps = lattice.get_symmetries()
# for rep in irreps:
#     print rep, irreps[rep]["k"], len(irreps[rep]["allowed_ops"])

# lattice.plot(coord_idx=False)
# lattice.plot_brillouinzone(loc="upper left")
# plt.show()

#############################
# Create sublattice structure
# DO NOT CHANGE ANYTHING HERE
#############################
sl1 = []
sl2 = []
sl3 = []
sl4 = []
span = np.array([[2.0, 0.0], [0.0, 2.0]])
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

    maxpg = ["C1", "C2", "C4", "C4^3"]
  
  
sublattice_structure = np.array([sl1,sl2,sl3,sl4])
print(sublattice_structure)

#############
#Define model
#############
# model = gm.GenericModel(lattice, maxpg=["C1", "C2", "C4s", "C4^3s"])

# model = gm.GenericModel(lattice, maxpg=maxpg, OBC=obc)
model = gm.GenericModel(lattice, OBC=obc)

model.set_nb_interaction(1, int_type='HOP', coupling='T')
model.set_nb_interaction(1, int_type='ISING', coupling='JZ')
model.set_nb_interaction(1, int_type='EXCHANGE', coupling='JXY')

# #HB with SCH
# model.set_interaction(basis_sites=[0,0,0], cellshifts=[[0,0],[1,0],[1,1]], int_type='GENERICSCH', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )
# model.set_interaction(basis_sites=[0,0,0], cellshifts=[[1,0],[1,1],[0,1]], int_type='GENERICSCH', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )
# model.set_interaction(basis_sites=[0,0,0], cellshifts=[[1,1],[0,1],[0,0]], int_type='GENERICSCH', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )
# model.set_interaction(basis_sites=[0,0,0], cellshifts=[[0,1],[0,0],[1,0]], int_type='GENERICSCH', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )

# model.plot()
# plt.show()


##################
# Create modelfile
##################
# model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".lat", 
#             sublattice_structure=sublattice_structure)
model.write(filename='square.' + str(lattice.n_sites) + "." + filename + ".lat")

# irreps = model.symmetries()
# irrepnames = irreps.keys()
# print "(\'" + "\' \'".join(irrepnames) + "\')"

