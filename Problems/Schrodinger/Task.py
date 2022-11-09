from Const import ab, dirac, m
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt

N = 300
a = ab
Xmin = -a
Xmax = a

omega = 10*dirac / (m*a**2)

def U(x):
    return m*(omega**2) * (x**2) / 2
    
def energy(n):
    E = dirac*omega*(n+0.5)
    b = np.sqrt(2*E/(m*omega**2))
    return E, -b, b

Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=4, mlt=0.2, mayRender=False, mayDrawEnergyLevels=False)
QuasiClassics.task(N, Xmin, Xmax, U, energy, ns = [0, 1, 2, 3], mlt=1, mayRender=False)
plt.legend()
plt.show()