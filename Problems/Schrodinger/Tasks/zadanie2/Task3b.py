import numpy as np
from Const import *
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.optimize

Z = 1
dirac = 1#(10**(-34))
m = 1#9.1 * (10**(-31))
ab = 1#0.5 * (10**(-10))
e = 1#1.6 * (10**(-19))

def A(a, b):
    return 1 / (np.sqrt(8 * ( 1/((2*a)**3 *(2*b)**3 ) + 1/(a+b)**6  )) * 4 *np.pi)

def T(a, b):
    return -dirac**2/m * (4*np.pi*A(a, b))**2 * ( -(a**2+b**2)/(16*a**3*b**3) - 8*a*b/(a+b)**6  )

def V(a, b):
    return -e**2 * Z * (4*np.pi*A(a, b))**2 * ( (a+b)/(8*a**3*b**3) + 8/(a+b)**5)

def U(a, b):
    return e**2 *(4*np.pi*A(a, b))**2 * (2.5 / (a+b)**5 + (a**2 + 3*a*b+b**2)/(8*a**2*b**2*(a+b)**3))

def N(a, b):
    return (4*np.pi*A(a, b))**2 * (1 / (8 * a**3 * b**3) + 8 / (a+b)**6)

def E(a, b):
    return (T(a, b) + V(a, b) + U(a, b))/N(a,b)

def Etrue():
    return -e**2/ab*(Z-5/16)**2

def alphaTrue():
    return (Z-5/16)/ab

print("true results: alpha={}, E={}".format(round(alphaTrue(), 2), round(Etrue(), 2)))

g = scipy.optimize.minimize(lambda arg: E(arg[0], arg[1]), (1, 1))
print("computed results: alpha={}, beta={}, E={}".format(round(g.x[0], 2), round(g.x[1], 2), round(g.fun, 2)))

b1 = 0.1
b2 = 3

bs = np.arange(b1, b2, (b2-b1)/100)

plt.plot(bs, [E(b, b) for b in bs])
plt.scatter([g.x[1]], [g.fun], color="red")
plt.show()
