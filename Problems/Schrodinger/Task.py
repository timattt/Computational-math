from Const import *
import numpy as np

N = 1000
a = 0.25*ab
Xmin = -1*a
Xmax = 1*a

def U(x):
    if np.abs(x) < a/2:
        return -100000 * e
    return 0
    
