#!/usr/bin/env python
import numpy as np
import quantipy.operators.heisenberg as heis

filename="bondmatrices.kitaevgamma.SU.2"
localdim=2

def get_kitaev_gamma_matrices(phi=0):
    """ phi = 0 -> Kitaev, phi = 0.5 -> Gamma """
    Sx, Sy, Sz = heis.spin_matrices(spin=2)
    x_matrix = np.sin(phi*np.pi)*(np.kron(Sy, Sz) + np.kron(Sz, Sy)) - \
               np.cos(phi*np.pi)*(np.kron(Sx, Sx))
    y_matrix = np.sin(phi*np.pi)*(np.kron(Sz, Sx) + np.kron(Sx, Sz)) - \
               np.cos(phi*np.pi)*(np.kron(Sy, Sy))
    z_matrix = np.sin(phi*np.pi)*(np.kron(Sx, Sy) + np.kron(Sy, Sx)) - \
               np.cos(phi*np.pi)*(np.kron(Sz, Sz))
    
    return x_matrix, y_matrix, z_matrix
    
print "write to file", filename

# Write Kitaev bonds
Kx, Ky, Kz, = get_kitaev_gamma_matrices(0)
fo = open(filename, "wb")
fo.write("[bondtype]=GENERICKITAEVX\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Kx:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")
fo.write("[bondtype]=GENERICKITAEVY\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Ky:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")
fo.write("[bondtype]=GENERICKITAEVZ\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Kz:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")


# Write Gamma bonds
Gx, Gy, Gz, = get_kitaev_gamma_matrices(.5)
fo.write("[bondtype]=GENERICGAMMAX\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Gx:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")
fo.write("[bondtype]=GENERICGAMMAY\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Gy:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")
fo.write("[bondtype]=GENERICGAMMAZ\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in Gz:
    fo.write(" ".join("({:1.3f},{:1.3f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

bondmatrix = heis.heisenberg_bond_sun(2, 1)
print bondmatrix
fo.write("[bondtype]=GENERICFOURSPIN\n")
fo.write("[sites]=4\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(bondmatrix, bondmatrix):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

# Write fourspin bond (diagonal)
fo.write("[bondtype]=GENERICFOURSPINDIAG\n")
fo.write("[sites]=4\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.diag(np.diag(np.kron(bondmatrix, bondmatrix))):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")



Sx, Sy, Sz = heis.spin_matrices(spin=2)


fo.write("[bondtype]=GENERICSXSX\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sx, Sx):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSXSY\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sx, Sy):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSXSZ\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sx, Sz):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSYSX\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sy, Sx):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSYSY\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sy, Sy):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSYSZ\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sy, Sz):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSZSX\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sz, Sx):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSZSY\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sz, Sy):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")

fo.write("[bondtype]=GENERICSZSZ\n")
fo.write("[sites]=2\n")
fo.write("[localdim]={0}\n".format(localdim))
for row in np.kron(Sz, Sz):
    fo.write(" ".join("({:1.5f},{:1.5f})".format(np.real(x),np.imag(x)) for x in row) + "\n")



fo.close()
