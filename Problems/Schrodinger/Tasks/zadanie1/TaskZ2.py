from Const import ab, dirac, m
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt
import math

N = 1000
a = ab
Xmin = 0
Xmax = 0.03*a

F = 1

def U(x):
    if x < 0:
        return 1000000000000000000000
    return F*x
    
def energy(n):
    E = math.pow(3*F*np.pi*dirac*(n+0.5)/(2*np.sqrt(2*m)), 2/3)
    return E, 0, E/F

Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=4, mlt=0.2, mayRender=False, mayDrawEnergyLevels=False)
QuasiClassics.task(N, Xmin, Xmax, U, energy, ns = [0, 1, 2, 3], mlt=1, mayRender=False, mayRound=False)
plt.legend()
plt.show()