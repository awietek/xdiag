#!/usr/bin/env python
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt

nsites = 20
J=1
Jd=1
seeds = range(1, 21, 1)

# TPQ
Ess = []
for seed in seeds:
    filename = "moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.seed.{}.txt".format(nsites, J, Jd, seed)
    data = np.loadtxt(filename)
    Ts = data[:,0]
    Es = data[:,2]
    plt.plot(Ts, Es, ":")
    Ess.append(Es)
Ess = np.array(Ess)
Es_mean = np.mean(Ess, axis=0)
Es_err = sem(Ess, axis=0)
plt.errorbar(Ts, Es_mean, Es_err, lw=3, c="tab:red", capsize=5, label="TPQ")
    
# Full ED
filename = "../checkerboard_fulled/moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.txt".format(nsites, J, Jd)
data = np.loadtxt(filename)
Ts = data[:,0]
Es = data[:,2]
plt.plot(Ts, Es, label="full ED", c="k", lw=2, zorder=100)
plt.xlabel(r"$T$")
plt.xlabel(r"$E$")
plt.legend()
plt.show()
