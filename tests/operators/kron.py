import numpy as np
from numpy import kron


sx = np.array([[0., 0.5], [0.5, 0.]])
sy = np.array([[0., -0.5j], [0.5j, 0.]])
sz = np.array([[0.5, 0.], [0.0, -0.5]])
jchimat = kron(sx, kron(sy, sz) - kron(sz, sy)) + \
          kron(sy, kron(sz, sx) - kron(sx, sz)) + \
          kron(sz, kron(sx, sy) - kron(sy, sx))

id = np.array([[1.0, 0.], [0.0, 1.0]])
id3 = kron(id, kron(id, id))


# p12 = (id3 + kron(kron(sx, sx) + kron(sy, sy) + kron(sz, sz), id) * 4)
# p23 = (id3 + kron(id, kron(sx, sx) + kron(sy, sy) + kron(sz, sz)) * 4)
# p123 = np.real(p23 @ p12 / 4)
# p132 = np.real(p12 @ p23 / 4)
# print(p123)
# print(p132)

# p12 = (id3 + kron(kron(sx, sx) + kron(sy, sy) + kron(sz, sz), id) * 4) / 4
# p23 = (id3 + kron(id, kron(sx, sx) + kron(sy, sy) + kron(sz, sz)) * 4) / 4
# p123 = np.real(p23 @ p12)
# p132 = np.real(p12 @ p23)
# # print(p123)
# # print(p132)

# # print(p123 - p132)
# print(np.imag(jchimat))
# # print((p123 - p132) + np.imag(jchimat) )


print(kron(sy, sz))
print(kron(sz, sy))
print(kron(sx, sz))
print(kron(sz, sx))
print(kron(sx, sy))
print(kron(sy, sx))
	       
