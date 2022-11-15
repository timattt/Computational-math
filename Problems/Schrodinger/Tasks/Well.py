from Const import ab, dirac, m, e
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt

N = 300
Xmin = -2*ab
Xmax = 2*ab

U0 = 6*(10**2)*e

def U(x):
    if abs(x) < ab:
        return -U0
    return 0

Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=5, mlt=1, mayDrawEnergyLevels=True)