from Const import ab, dirac, m
import numpy as np
import Schrodinger
import QuasiClassics
import matplotlib.pyplot as plt


def error(N):
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
    
    Es, psis, total = Schrodinger.task(N, Xmin, Xmax, U, maxDiscrete=4, mlt=0.2, mayRender=False, mayDrawEnergyLevels=False, mayPrint=False)
    plt.clf()
    err = 0
    for i in range(3):
        tE, tmp1, tmp2 =  energy(i)
        err += abs(Es[i] - tE)
        
    return err
    
    
ns = []
errs = []
    
for N in range(10, 20):
    ns.append(np.log(N))
    errs.append(np.log(error(N)))
    
coefs = np.polyfit(ns, errs, 1)
print("order={}".format(-int(coefs[0])))
plt.plot(ns, [coefs[0] * n + coefs[1] for n in ns], label = "solution")
plt.plot(ns, errs, label = "linear aprox")
plt.xlabel("log(N)")
plt.ylabel("log(err)")
plt.legend()
plt.grid()
plt.show()