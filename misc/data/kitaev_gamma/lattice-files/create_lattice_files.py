#!/usr/bin/env python
import quantipy.lattice as latt
import quantipy.models.genericmodel as gm
import quantipy.quicked.quicked as qed
import quantipy.operators.heisenberg as heis
import quantipy.spectra.spectra as spec
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.sparse.linalg
import sys

filename = "HeisenbergKitaevGamma"

lattice = latt.Honeycomb(simulation_torus_matrix=[[2, 0], [0, 3]])

# lattice = latt.Honeycomb(simulation_torus_matrix=[[1, 4], [5, -1]])   # 42 site
# lattice = latt.Honeycomb(simulation_torus_matrix=[[4, 0], [0, 4]])   # 32 site
# lattice = latt.Honeycomb(simulation_torus_matrix=[[2, 2], [2, -4]])   # 24 site
# lattice = latt.Honeycomb(simulation_torus_matrix=[[3, 0], [0, 3]])   # 18 site
# lattice = latt.Honeycomb(simulation_torus_matrix=[[2, 0], [0, 2]])    # 8 site

model = gm.GenericModel(lattice)

#############################
# Create sublattice structure
# DO NOT CHANGE ANYTHING HERE
#############################
a1=lattice.a1
a2=lattice.a2
sl1 = [] 
sl2 = [] 
sl3 = [] 
sl4 = [] 
for idx, coord in enumerate(lattice.coordinates):
    A = np.array([a1, 2*a2]).T
    c = np.linalg.solve(A,coord)
    case = np.array(np.round(c*2),dtype=np.int) %2
    print(idx, coord, case)
    if (case == np.array([0,0])).all():
        sl1.append(idx)
    elif (case == np.array([1,0])).all():
        sl2.append(idx)
    elif (case == np.array([0,1])).all():
        sl3.append(idx)
    elif (case == np.array([1,1])).all():
        sl4.append(idx)

    else:
        print("somethings wrong")

sublattice_structure = np.array([sl1, sl2, sl3, sl4])
print(sl1)
print(sl2)
print(sl3)
print(sl4)
if lattice.n_sites == 42:
    sublattice_structure = np.array([sl1 + sl3, sl2 + sl4])
print(sublattice_structure)
model.set_nb_interaction(0, int_type='GENERICSZ', coupling='HZ')

# # Define Heisenberg model
model.set_nb_interaction(1, int_type='HB', coupling='J')


# Define Gamma model
model.set_interaction(np.array([[0,0,0],[0,0,1]]),
                      int_type='GENERICKITAEVX', coupling='K',
                      directed=False)
model.set_interaction(np.array([[1,-1,1],[0,0,0]]),
                      int_type='GENERICKITAEVY', coupling='K',
                      directed=False)
model.set_interaction(np.array([[1,0,1],[0,0,0]]),
                      int_type='GENERICKITAEVZ', coupling='K',
                      directed=False)
model.set_interaction(np.array([[0,0,0],[0,0,1]]),
                      int_type='GENERICGAMMAX', coupling='G',
                      directed=False)
model.set_interaction(np.array([[1,-1,1],[0,0,0]]),
                      int_type='GENERICGAMMAY', coupling='G',
                      directed=False)
model.set_interaction(np.array([[1,0,1],[0,0,0]]),
                      int_type='GENERICGAMMAZ', coupling='G',
                      directed=False)

model.plot(coord_idx=True)
model.plot_brillouinzone()
plt.show()

model.write(filename='honeycomb.'+str(lattice.n_sites) + "." + filename + ".fsl.lat", 
            sublattice_structure= sublattice_structure)
# for phi in [0, 0.25, 0.3, 0.35 ,0.4, 0.45 ,0.5 ]:
#     print phi, -np.cos(phi*np.pi), np.sin(phi*np.pi)

# hb_matrix = heis.heisenberg_bond(spin=2)
# type_matrices = {'HB': hb_matrix}
# couplings = {'J1': 1.}
# bonds = qed.get_bonds_from_model(model, type_matrices, couplings)

# # Get bondmatrices and coupling strengths for whole moded
# Sx, Sy, Sz = heis.spin_matrices(spin=2)
# x_matrix = np.sin(phi*np.pi)*(np.kron(Sy, Sz) + np.kron(Sz, Sy)) - \
#            np.cos(phi*np.pi)*(np.kron(Sx, Sx))
# y_matrix = np.sin(phi*np.pi)*(np.kron(Sz, Sx) + np.kron(Sx, Sz)) - \
#            np.cos(phi*np.pi)*(np.kron(Sy, Sy))
# z_matrix = np.sin(phi*np.pi)*(np.kron(Sx, Sy) + np.kron(Sy, Sx)) - \
#            np.cos(phi*np.pi)*(np.kron(Sz, Sz))
# print x_matrix
# # x_matrix = np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)
# # y_matrix = np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)
# # z_matrix = np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)
# np.set_printoptions(precision=20)
# type_matrices = {'Gammax': x_matrix, 'Gammay': y_matrix, 'Gammaz': z_matrix}
# couplings = {'Jx': 1., 'Jy': 1., 'Jz': 1.}
