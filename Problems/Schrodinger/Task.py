from Const import *
import numpy as np

N = 1000
a = 0.25*ab
Xmin = 1*ab
Xmax = 25*ab

A = 5*ab
U0 = 1*e

def U(x):
    return U0 * (x/A - A/x)**2
    
