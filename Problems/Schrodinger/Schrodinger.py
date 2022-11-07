import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
import math
from Const import *

import Task

# init values
h = (Task.Xmax - Task.Xmin) / Task.N
xs = np.array([Task.Xmin + i *h for i in range(Task.N)])
U = np.array([Task.U(x) for x in xs])
Umin = np.min(U)
Umax = np.max(U)
solid = -(dirac**2) / (2*m)
beta = solid / (h**2)
alpha = np.array([U[i] - 2*solid / (h**2) for i in range(Task.N)])
A = np.zeros((Task.N, Task.N))

# init matrix
for i in range(Task.N):
    A[i, i] = alpha[i]
    
for i in range(Task.N-1):
    A[i, i+1] = A[i+1, i] = beta

# solve spectral problem
w, v = np.linalg.eigh(A)

print("Umin = {}, Umax = {}".format(Umin, Umax))

# normalize
psis = np.array([v[i] / np.sqrt(np.dot(v[i], v[i])*h) for i in range(Task.N)])
psiMax = np.max(np.abs(psis))
totalDiscrete = np.sum([1 if w[i] < Umax and w[i] > Umin else 0 for i in range(Task.N)])

xs = xs

for i in range(Task.N):
        psi_ = psis[i]
        E = w[i]
        if E < Umax and E > Umin:
            psi_ *= (Umax-Umin)/(np.max(psi_)-np.min(psi_))/e/totalDiscrete
            psi_ += E/e
            plt.plot(xs/ab, psi_/1000000, label="E={}".format(E/e))

plt.plot(xs/ab, U/e/1000000)            
plt.xlabel("радиус бора")
plt.ylabel("МЭВ")
print("Total discrete levels: {}".format(totalDiscrete))
plt.legend()
plt.show()