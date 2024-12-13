#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import quantipy.lattice as latt
import quantipy.models.genericmodel as gm
import scipy
import scipy.linalg as linalg
from matplotlib import rc
rc('font',**{'family':'Times New Roman', 'size': 13})
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}')


##################################################
# Script to create model files for Kagome lattices
# for the EDFSL engine.
# prefix of filename is "kagome.$n_sites"
# postfix is ".fsl.lat"
##################################################

filename = "J1J2.sublattices" 

a1 = np.array([1., 0.])
a2 = np.array([.5, np.sqrt(3)/2])

# 12 sites 
lattice = latt.Triangular(simulation_torus=np.array([2*a1 + 2*a2, 2*a1 - 4*a2]))
filename = "J1J2" 

# # 18 sites 
# lattice = latt.Triangular(simulation_torus=np.array([1*a1 + 4*a2, 4*a1 - 2*a2]))
# filename = "10418.J1J2.sublattices" 

# # 21 sites
# lattice = latt.Triangular(simulation_torus=np.array([1*a1 + 4*a2, 5*a1 - 1*a2]))
# filename = "10421.J1J2.sublattices" 

# # 24 sites 
# lattice = latt.Triangular(simulation_torus=np.array([1*a1 + 4*a2, 5*a1 - 4*a2]))
# filename = "10424.J1J2.sublattices" 

# ########################
# # Print and Plot lattice
# ########################

# # for idx, coord in enumerate(lattice.coordinates):
# #     print idx, coord
# # perms, irreps = lattice.get_symmetries()
# # for rep in irreps:
# #     print rep, irreps[rep]["k"], len(irreps[rep]["allowed_ops"])


lattice.plot(coord_idx=True)
plt.show()
# plt.savefig('triangular.'+str(lattice.n_sites) + "." + filename + ".tsl.pdf")
# lattice.plot_brillouinzone()
# plt.savefig('triangular.'+str(lattice.n_sites) + "." + filename + ".tsl.bz.pdf")



#############################
# Create sublattice structure
# DO NOT CHANGE ANYTHING HERE
#############################
sl1 = [] 
sl2 = [] 
sl3 = [] 
for idx, coord in enumerate(lattice.coordinates):
    A = np.array([a1+a2, 3*a2]).T
    c = linalg.solve(A,coord)
    case = np.array(np.round(c*3),dtype=np.int64) %3
    # print idx, coord, case
    if (case == np.array([0,0])).all():
        sl1.append(idx)
    elif (case == np.array([0,1])).all():
        sl2.append(idx)
    elif (case == np.array([0,2])).all():
        sl3.append(idx)
    else:
        # print "somethings wrong"
        pass
sublattice_structure = np.array([sl1, sl2, sl3])
# print sublattice_structure



#############
#Define model
#############
# model = gm.GenericModel(lattice, maxpg=["C1"])
model = gm.GenericModel(lattice)

# model.set_nb_interaction(1)
# model.set_interaction(np.array([[1,1,0],[0,1,1],[1,0,2]]), int_type='ScalarChirality', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )
# model.set_interaction(np.array([[0,0,0],[0,0,1],[0,0,2]]), int_type='ScalarChirality', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )


# #TFI # model.set_nb_interaction(1, int_type='SzSz', coupling='J')
# model.set_nb_interaction(0, int_type='HX', coupling='H')

#J1J2
model.set_nb_interaction(1, int_type='SdotS', coupling='J1')
model.set_nb_interaction(2, int_type='SdotS', coupling='J2')

# # ScalarChirality
# model.set_interaction(np.array([[0,0,0],[1,0,0],[0,1,0]]), int_type='ScalarChirality', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )
# model.set_interaction(np.array([[0,0,0],[0,1,0],[-1,1,0]]), int_type='ScalarChirality', coupling='Jchi', directed=True, allowed_perms=np.array([[0,1,2],[1,2,0],[2,0,1]]) )

# model.plot_interactions()
# plt.show()
##################
# Create modelfile
##################

# model.write(filename='triangular.'+str(lattice.n_sites) + "." + filename + ".tsl.lat", 
#                        sublattice_structure= sublattice_structure)

# print "('" + "' '".join(model._do_symmetries().keys()) +"')"
model.write(filename='triangular.'+str(lattice.n_sites) + "." + filename + ".lat")

