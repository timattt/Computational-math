from Const import ab, dirac, m
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt
import StationaryPerturbationTheory

N = 100
a = 10*ab
Xmin = -a
Xmax = a

omega = 10*dirac / (m*a**2)
A = 10000000000
B = 10000000000000000000#19

def V(x):
    return A * x**3 + B * x**4

def U(x):
    return m*(omega**2) * (x**2) / 2 + V(x)

StationaryPerturbationTheory.task(N, Xmin, Xmax, U, V, maxDiscrete = 4, mlt=1, mayRound = False)