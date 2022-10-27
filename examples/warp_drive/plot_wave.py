#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

time = 60

warps = [1, 2, 4]
for wf in warps:
    data = np.loadtxt("output.warp.factor.{:.1f}.txt".format(wf))[time, :]
    plt.plot(data, label=r"$\Omega={}$".format(wf))

plt.legend()
plt.show()
    
