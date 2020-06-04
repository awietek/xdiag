#!/usr/bin/env python
import numpy as np

N=6
bonds = [(i, (i+1)%N) for i in range(N)]

mat = np.zeros((N,N))

for bond in bonds:
    s1 = bond[0]
    s2 = bond[1]
    mat[s1, s2] = 1.
    mat[s2, s1] = 1.

print(mat)
print(np.linalg.eigvalsh(mat))


