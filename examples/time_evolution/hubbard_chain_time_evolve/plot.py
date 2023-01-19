#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("output.txt")
im = plt.imshow(data, aspect='auto', cmap="cividis",
                extent=(0, 10, 0, 10), origin = "lower")
plt.colorbar(im)
plt.title(r"Density evolution of a Hubbard doublon, $\langle n_i(t) \rangle$")
plt.xlabel("i")
plt.ylabel("t")
plt.savefig("hubbard_doublon.png")
plt.show()


