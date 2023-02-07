#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

n_sites = 20
J=1
Jd=1
seeds = range(1, 3, 1)

# TPQ
for seed in seeds:
    filename = "moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.seed.{}.txt".format(n_sites, J, Jd, seed)
    data = np.loadtxt(filename)
    Ts = data[:,0]
    Es = data[:,2]
    plt.plot(Ts, Es)

# Full ED
filename = "../checkerboard_fulled/moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.txt".format(n_sites, J, Jd)
data = np.loadtxt(filename)
Ts = data[:,0]
Es = data[:,2]
plt.plot(Ts, Es, label="full ED", c="k")
plt.xlabel(r"$T$")
plt.xlabel(r"$E$")
plt.legend()
plt.show()
