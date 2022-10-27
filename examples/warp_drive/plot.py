#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

warp_factor = 0.0

data = np.loadtxt("output.warp.factor.{:.1f}.txt".format(warp_factor))
im = plt.imshow(data, aspect='auto', cmap="seismic")
plt.colorbar(im)
plt.title(r"$\Omega={}$".format(warp_factor))
plt.savefig("output.warp.factor.{:.1f}.pdf".format(warp_factor))
plt.show()


