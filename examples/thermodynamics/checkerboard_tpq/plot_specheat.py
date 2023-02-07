#!/usr/bin/env python
import numpy as np
from scipy.stats import sem
import matplotlib.pyplot as plt

n_sites = 20
J=1
Jd=1
seeds = range(1, 21, 1)

# TPQ
Ess = []
Qss = []
for seed in seeds:
    filename = "moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.seed.{}.txt".format(n_sites, J, Jd, seed)
    data = np.loadtxt(filename)
    Ts = data[:,0]
    Es = data[:,2]
    Qs = data[:,3]
    plt.plot(Ts, (1 / Ts)**2 * (Qs - Es**2), ":")
    Ess.append(Es)
    Qss.append(Qs)

Ess = np.array(Ess)
Qss = np.array(Qss)
        
def jackknife(data):
    """ Resample to Jackknife averages """
    data_resampled = np.zeros_like(data)
    for idx in range(data.shape[0]):
        data_resampled[idx, :] = np.mean(np.delete(data, idx, axis=0), axis=0)
    return data_resampled

Ess = jackknife(Ess)
Qss = jackknife(Qss)
Css = (Qss - Ess**2) / (Ts**2)
Cs_mean = np.mean(Css, axis=0)
Cs_error = (Css.shape[0]-1) * sem(Css, axis=0)
plt.errorbar(Ts, Cs_mean, Cs_error, lw=3, c="tab:red", capsize=5, label="TPQ")
    
# Full ED
filename = "../checkerboard_fulled/moments/moments.checkerboard.{}.J.{:.2f}.Jd.{:.2f}.txt".format(n_sites, J, Jd)
data = np.loadtxt(filename)
Ts = data[:,0]
Cs = data[:,4]
plt.plot(Ts, Cs, label="full ED", c="k", lw=2, zorder=100)
plt.xlabel(r"$T$")
plt.xlabel(r"$C$")
plt.legend()
plt.show()
