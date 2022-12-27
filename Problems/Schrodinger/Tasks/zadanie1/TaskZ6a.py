from Const import ab, dirac, m
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt
import StationaryPerturbationTheory

N = 100
a = ab
Xmin = -a
Xmax = a

omega = 10*dirac / (m*a**2)
alpha = 0.000001

def V(x):
    return alpha * x

def U(x):
    return m*(omega**2) * (x**2) / 2 + V(x)

StationaryPerturbationTheory.task(N, Xmin, Xmax, U, V, maxDiscrete = 4, mlt=1)