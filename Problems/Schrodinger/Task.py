from Const import *
import numpy as np
import Schrodinger

N = 1000
Xmin = -16000000*ab
Xmax = 16000000*ab

omega = 100*dirac / (m*np.sqrt(ab))

def U(x):
    return m*omega**2 * x**2 / 2
    
Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=4)
